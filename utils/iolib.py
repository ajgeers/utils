"""Module for reading, writing, compressing and converting files.

Please note, some of the functions in this module were created and tested using
VTK 5. VTK 6 introduced a number of backwards-incompatible changes, including
replacing 'SetInput()' with 'SetInputData()' and 'SetInputConnection'.

"""
import glob
import os
import csv
import vtk
import gzip
import StringIO


def compress(path='test.vtp'):
    """Compress file with gzip."""
    with open(path, 'rb') as ifile:
        with gzip.open(path + '.gz', 'wb') as ofile:
            ofile.writelines(ifile)


def decompress(path='test.vtp.gz'):
    """Decompress file with gzip."""
    with gzip.open(path, 'rb') as ifile:
        with open(path[:-3], 'wb') as ofile:
            ofile.write(ifile.read())


def csv_to_list(path):
    """Convert CSV-file to a nested list of strings."""
    with open(path, 'rb') as f:
        reader = csv.reader(f)
        return list(reader)


def csv_to_dict(path):
    """Create nested dictionary from csv file. Workaround for when pandas is
    unavailable and you want to select 2D array elements with row and column
    names rather than integers.

    * First row is used for column names
    * First column is used for row names.
    * Access data from dictionary x using x['rowname']['columnname']
    * Extract all row names with x.keys()
    * Extract all column names with x.values()[0].keys()

    Note: Expects '\n' as newline character.

    """
    x = {}
    with open(path, 'rb') as f:
        header = f.next().strip().split(',')[1:]
        for line in f:
            row = line.strip().split(',')
            x[row[0]] = dict(
                (header[i], v) for i, v in enumerate(row[1:]))
    return x


def listdir(path, match='*', dirname=False, extension=False):
    """List all files and folders in specified directory.

    Args:
        path: Path to directory.
        match: Specify file name pattern according to rules used by Unix
            shell. For instance, 'match=*.pdf' gives you a list of names of all
            the pdf-files in 'path'.
        dirname (bool): Include whole path name.
        extension (bool): Include file extension.

    """
    items = glob.glob(os.path.join(path, match))
    if not dirname:
        items = [os.path.basename(item) for item in items]
    if not extension:
        items = [os.path.splitext(item)[0] for item in items]
    return items


def readvti(path):
    """Read VTI-file, i.e. image in VTK XML format."""
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutput()


def readvtk(path, datatype='polydata'):
    """Read VTK-file.

    Args:
        path: Path to file.
        type: 'imagedata', 'polydata', 'unstructeredgrid'

    """
    if datatype=='imagedata':
        reader = vtk.vtkImageDataReader()
    elif datatype=='polydata':
        reader = vtk.vtkPolyDataReader()
    elif datatype=='unstructeredgrid':
        reader = vtk.vtkUnstructuredGridReader()
    else:
        print 'Invalid datatype'
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutput()


def readvtp(path, dataarrays=True):
    """Read VTP-file, i.e. polydata in VTK XML format.

    Args:
        dataarrays (bool): Include point and cell data.

    """
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    if dataarrays == False:
        for i in range(reader.GetNumberOfPointArrays()):
            arrayname = reader.GetPointArrayName(i)
            reader.SetPointArrayStatus(arrayname, 0)
        for i in range(reader.GetNumberOfCellArrays()):
            arrayname = reader.GetCellArrayName(i)
            reader.SetPointArrayStatus(arrayname, 0)
        reader.Update()
    return reader.GetOutput()


def readvtu(path):
    """Read VTU-file, i.e. unstructured grid in VTK XML format."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutput()


def replacestring(lines, tag, value):
    """Replace string in list of strings.

    Args:
        lines: List of strings.
        tag: String to replace.
        value: String with which to replace 'tag'.

    """
    output = []
    for line in lines:
        line = line.replace(tag, value)
        output.append(line)
    return output


def writepoints(points, filename):
    """Write points as VTP-file."""
    polydata = vtk.vtkPolyData()
    cellarray = vtk.vtkCellArray()

    for i in range(points.GetNumberOfPoints()):
        cellarray.InsertNextCell(1)
        cellarray.InsertCellPoint(i)

    polydata.SetPoints(points)
    polydata.SetVerts(cellarray)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInput(polydata)
    writer.Write()


def writevti(image, path):
    """Write VTI-files, i.e. images in VTK XML format."""
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetInput(image)
    writer.SetFileName(path)
    writer.Write()


def writevtp(polydata, path):
    """Write VTP-files, i.e. polydata in VTK XML format."""
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInput(polydata)
    writer.SetFileName(path)
    writer.Write()

def writevtu(grid, path):
    """Write VTU-files, i.e. unstructured grids in VTK XML format."""
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInput(grid)
    writer.SetFileName(path)
    writer.Write()


#-------------------------------------------------------------------------------
# CFX
#-------------------------------------------------------------------------------

def cfx2vtp(inputfile, outputfile, surface=True, ascii=False):
    """Convert polydata exported from CFX-Post to VTP.

    Args:
        surface (bool): Convert surface or line polydata.
        ascii (bool): Return VTP file in ASCII format.

    Export surface in CFX-Post with following options:
    * file extension: csv
    * export geometry information: line and face connectivity
    * (optional) select variable(s)
    * vector display: scalar
    * separator: comma space
    * include header

    """
    f = open(inputfile, 'rb')

    # derive data size from csv file
    if surface:
        for i, line in enumerate(f):
            if line.strip() == '[Data]':
                datalinenumber = i
            if line.strip() == '[Faces]':
                faceslinenumber = i
            lastlinenumber = i
        numberofnodes = faceslinenumber - datalinenumber - 3
        numberofelements = lastlinenumber - faceslinenumber - 1
    else:
        for i, line in enumerate(f):
            if line.strip() == '[Data]':
                datalinenumber = i
            if line.strip() == '[Lines]':
                lineslinenumber = i
        numberofnodes = lineslinenumber - datalinenumber - 3

    # obtain list of variables names
    f.seek(0)
    for i in range(datalinenumber + 2):
        arrayline = f.readline()
    arraynames = arrayline.strip().split(', ')
    arraynames[0:3] = []

    # define polydata
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    points.SetNumberOfPoints(numberofnodes)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells) if surface else polydata.SetLines(cells)
    for arrayname in arraynames:
        array = vtk.vtkDoubleArray()
        array.SetName(arrayname)
        array.SetNumberOfTuples(numberofnodes)
        polydata.GetPointData().AddArray(array)

    # parse through the rest of the file using the csv module
    reader = csv.reader(f)

    # assign x,y,z coordinates and variable values to points
    for i in range(numberofnodes):
        dataline = reader.next()
        point = [float(dataline[0]), float(dataline[1]), float(dataline[2])]
        points.SetPoint(i, point)
        for j in range(len(arraynames)):
            dataarray = polydata.GetPointData().GetArray(arraynames[j])
            dataarray.SetComponent(i, 0, float(dataline[j + 3]))

    # skip element '[Faces]' (or '[Lines]') in csv-file
    reader.next()
    reader.next()

    if surface:
        # obtain and set connectivity
        cellids = vtk.vtkIdList()
        for i in range(numberofelements):
            facesline = reader.next()
            cellids.Initialize()
            for j in range(len(facesline)):
                cellids.InsertNextId(int(facesline[j]))
            cells.InsertNextCell(cellids)
    else:
        # obtain connectivity
        connectivitylist = []
        for row in reader:
            row = [int(item) for item in row]
            connectivitylist.append(row)
        connectivitylist = filter(None, connectivitylist)

        # rearrange connectivity
        linecounter = 0
        for i in range(len(connectivitylist)):
            if i == 0:
                connectivity = [connectivitylist[i]]
            elif connectivitylist[i][0] == connectivitylist[i - 1][1]:
                connectivity[linecounter].append(connectivitylist[i][1])
            else:
                connectivity.append([])
                linecounter += 1
                connectivity[linecounter].append(connectivitylist[i][0])
                connectivity[linecounter].append(connectivitylist[i][1])

        # set connectivity
        cellids = vtk.vtkIdList()
        for i in range(len(connectivity)):
            cellids.Initialize()
            for j in range(len(connectivity[i])):
                cellids.InsertNextId(int(connectivity[i][j]))
            cells.InsertNextCell(cellids)

    f.close()

    # write vtk polydata
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInput(polydata)
    if ascii: writer.SetDataModeToAscii()
    writer.SetFileName(outputfile)
    writer.Write()


def vtp2cfx(inputfile, outputfile, surface=True):
    """Convert VTP polydata to format that can be imported into CFX-Post.

    Args:
        surface (bool): Convert surface or line polydata.

    """
    # read vtp file
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(inputfile)
    reader.Update()
    polydata = reader.GetOutput()

    # read names of data arrays
    arraynames = []
    dataarrays = polydata.GetPointData()
    numberofdataarrays = dataarrays.GetNumberOfArrays()
    for i in range(numberofdataarrays):
        array = dataarrays.GetArray(i)
        arrayname = array.GetName()
        arraynames.append(arrayname)

    # append names of data arrays to header and write header
    f = open(outputfile, 'wb')
    header = "\n[Name]\nSEGMENT\n\n[Data]\nX [ m ], Y [ m ], Z [ m ]"
    for i in range(numberofdataarrays):
        header += ", " + arraynames[i]
    header += "\n"
    f.write(header)

    # write values of x,y,z and data arrays row by row
    for i in range(polydata.GetNumberOfPoints()):
        point = polydata.GetPoint(i)
        line = str(point[0]) + ', ' + str(point[1]) + ', ' + str(point[2])
        for arrayname in arraynames:
            array = dataarrays.GetArray(arrayname)
            line += ', ' + str(array.GetComponent(i, 0))
        line += '\n'
        f.write(line)

    # write list of connectivity
    if surface:
        line = '\n[Faces]\n'
        f.write(line)
        for i in range(polydata.GetNumberOfCells()):
            cellpointids = polydata.GetCell(i).GetPointIds()
            line = ''
            for j in range(cellpointids.GetNumberOfIds()):
                if (j > 0):
                    line += ', '
                line += str(cellpointids.GetId(j))
            line += '\n'
            f.write(line)
    else:
        line = '\n[Lines]\n'
        f.write(line)
        for i in range(polydata.GetNumberOfCells()):
            cellpointids = polydata.GetCell(i).GetPointIds()
            line = ''
            for j in range(cellpointids.GetNumberOfIds() - 1):
                line += (str(cellpointids.GetId(j)) + ', ' +
                         str(cellpointids.GetId(j + 1)) + '\n')
            f.write(line)

    # add blank line to mimic exact same file structure as CFX-generated
    # csv-file
    line = '\n'
    f.write(line)

    f.close()
