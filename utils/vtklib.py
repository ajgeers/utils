"""Module for operating on VTK data objects.

For more details on VTK, go to www.vtk.org. For more details specifically on VTK
classes, go to http://www.vtk.org/doc/release/5.10/html/classes.html.

Please note, the functions in this module were created and tested using VTK
5.10. VTK 6 introduced a number of backwards-incompatible changes, including
replacing 'SetInput()' with 'SetInputData()' and 'SetInputConnection'.

"""
import math
import vtk


def addtoarray(polydata, inputarrayname, outputarrayname, value):
    """Add value to inputarray of polydata and store in outputarray."""
    calc = vtk.vtkArrayCalculator()
    calc.SetInput(polydata)
    calc.AddScalarArrayName(inputarrayname, 0)
    calc.SetFunction('{0} + {1}'.format(inputarrayname, value))
    calc.SetResultArrayName(outputarrayname)
    calc.Update()
    return calc.GetOutput()


def appendfilter(polydata1, polydata2):
    """Append two polydatas."""
    appender = vtk.vtkAppendPolyData()
    appender.AddInput(polydata1)
    appender.AddInput(polydata2)
    appender.Update()
    return appender.GetOutput()


def areaweightedmean(surface, arrayname, scalarabsolute=False):
    """Compute the area-weighted mean of a scalar or vector field.

    Args:
        surface: Surface mesh.
        arrayname: Name scalar or vector field.
        scalarabsolute (bool): Take absolute value of scalar field before
            averaging. Has no effect if array is a vector field.

    Returns:
        Area-averaged scalar or vector.

    Note:
        Pointdata is first converted to celldata, because it is easier to obtain
        the surface area represented by a cell than by a point. Then, we loop
        through all cells and compute the average value of the scalar or vector
        field weighted by the surface area of the cells.

    """
    surface = pointdatatocelldata(surface)  # convert pointdata to celldata
    array = surface.GetCellData().GetArray(arrayname)
    isscalar = True if array.GetNumberOfComponents() == 1 else False
    sumarea = 0
    if isscalar:  # array is a scalar field
        sumscalar = 0
        for i in range(surface.GetNumberOfCells()):
            area = surface.GetCell(i).ComputeArea()
            if scalarabsolute:
                scalar = abs(array.GetValue(i))
            else:
                scalar = array.GetValue(i)
            sumarea += area
            sumscalar += area * scalar
        mean = sumscalar / sumarea
    else:  # array is a vector field
        sumvector = [0.0, 0.0, 0.0]
        for i in range(surface.GetNumberOfCells()):
            area = surface.GetCell(i).ComputeArea()
            vector = [array.GetComponent(i, 0),
                      array.GetComponent(i, 1),
                      array.GetComponent(i, 2)]
            sumarea += area
            sumvector = [sumvector[0] + area * vector[0],
                         sumvector[1] + area * vector[1],
                         sumvector[2] + area * vector[2]]
        mean = [sumvector[0] / sumarea,
                sumvector[1] / sumarea,
                sumvector[2] / sumarea]
    return mean


def booldifference(surface1, surface2):
    """Create surface that is the difference between the two input surfaces."""
    boolfilter = vtk.vtkBooleanOperationPolyDataFilter()
    boolfilter.SetOperationToDifference()
    boolfilter.SetInput(0, surface1)
    boolfilter.SetInput(1, surface2)
    boolfilter.Update()
    return boolfilter.GetOutput()


def celldatatopointdata(polydata):
    """Convert celldata to pointdata."""
    converter = vtk.vtkCellDataToPointData()
    converter.SetInput(polydata)
    converter.Update()
    return converter.GetOutput()


def cellnormals(surface):
    """Add cell normals."""
    normals = vtk.vtkPolyDataNormals()
    normals.SetInput(surface)
    normals.ComputePointNormalsOff()
    normals.ComputeCellNormalsOn()
    normals.Update()
    return normals.GetOutput()


def centerofmass(surface):
    """Compute the center of mass of a surface mesh.

    The center of mass of the surface mesh is computed by area-weighted
    averaging the center of mass (i.e. centroid) of all triangles in the
    mesh.

    """
    surface = triangulate(surface)  # ensure elements are triangles
    sumarea = 0
    sumpoint = [0.0, 0.0, 0.0]
    for i in range(surface.GetNumberOfCells()):
        currentcell = surface.GetCell(i)
        area = currentcell.ComputeArea()
        point0 = [currentcell.GetPoints().GetPoint(0)[0],
                  currentcell.GetPoints().GetPoint(0)[1],
                  currentcell.GetPoints().GetPoint(0)[2]]
        point1 = [currentcell.GetPoints().GetPoint(1)[0],
                  currentcell.GetPoints().GetPoint(1)[1],
                  currentcell.GetPoints().GetPoint(1)[2]]
        point2 = [currentcell.GetPoints().GetPoint(2)[0],
                  currentcell.GetPoints().GetPoint(2)[1],
                  currentcell.GetPoints().GetPoint(2)[2]]
        point = [(point0[0] + point1[0] + point2[0]) / 3.0,
                 (point0[1] + point1[1] + point2[1]) / 3.0,
                 (point0[2] + point1[2] + point2[2]) / 3.0]
        sumarea += area
        sumpoint = [sumpoint[0] + area * point[0],
                    sumpoint[1] + area * point[1],
                    sumpoint[2] + area * point[2]]
    meanpoint = [sumpoint[0] / sumarea,
                 sumpoint[1] / sumarea,
                 sumpoint[2] / sumarea]
    return meanpoint


def cleanpolydata(polydata):
    """Apply VTK mesh cleaning filter to polydata."""
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInput(polydata)
    cleaner.Update()
    return cleaner.GetOutput()


def countregions(polydata):
    """Count number of disconnected regions in polydata."""
    connect = vtk.vtkPolyDataConnectivityFilter()
    connect.SetInput(polydata)
    connect.Update()
    return connect.GetNumberOfExtractedRegions()


def delaunay2d(points):
    """Construct a 2D Delaunay triangulation from a set of points."""
    delaunay = vtk.vtkDelaunay2D()
    delaunay.SetInput(points)
    delaunay.Update()
    return delaunay.GetOutput()


def pointdistance(point1, point2):
    """Compute Euclidean distance between two points."""
    return math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point1, point2))


def extractboundaryedge(surface, feature_edges=False):
    """Extract boundary edges of a surface mesh."""
    edge = vtk.vtkFeatureEdges()
    edge.SetInput(surface)
    if not feature_edges:
        edge.FeatureEdgesOff()
    edge.Update()
    return edge.GetOutput()


def extractcells(polydata, idlist=[0, 1, 2]):
    """Extract cells from polydata whose cellid is in idlist."""
    cellids = vtk.vtkIdList()  # specify cellids
    cellids.Initialize()
    for i in idlist:
        cellids.InsertNextId(i)
    extract = vtk.vtkExtractCells()  # extract cells with specified cellids
    extract.SetInput(polydata)
    extract.AddCellList(cellids)
    extraction = extract.GetOutput()
    geometry = vtk.vtkGeometryFilter()  # unstructured grid to polydata
    geometry.SetInput(extraction)
    return geometry.GetOutput()


def extractclosestpointregion(polydata, point=[0, 0, 0]):
    """Extract region closest to specified point."""
    connect = vtk.vtkPolyDataConnectivityFilter()
    connect.SetInput(polydata)
    connect.SetExtractionModeToClosestPointRegion()
    connect.SetClosestPoint(point)
    connect.Update()
    return connect.GetOutput()


def extractedges(surface):
    """Extract edges of a surface mesh, i.e. the edges of all mesh elements."""
    extractor = vtk.vtkExtractEdges()
    extractor.SetInput(surface)
    extractor.Update()
    return extractor.GetOutput()


def extractlargestregion(polydata):
    """Extract largest of several disconnected regions."""
    connect = vtk.vtkPolyDataConnectivityFilter()
    connect.SetInput(polydata)
    connect.SetExtractionModeToLargestRegion()
    connect.Update()
    return connect.GetOutput()


def extractsurface(unstructuredgrid):
    """Extract surface mesh from an vtkUnstructuredGrid object."""
    surfer = vtk.vtkDataSetSurfaceFilter()
    surfer.SetInput(unstructuredgrid)
    surfer.Update()
    return surfer.GetOutput()


def fillholes(surface, holesize=1000):
    """Fill holes in surface. Use holesize to specify the max 'radius' of the
    holes to be filled."""
    filler = vtk.vtkFillHolesFilter()
    filler.SetInput(surface)
    filler.SetHoleSize(holesize)
    filler.Update()
    return filler.GetOutput()


def findclosestpoint(polydata, refpoint):
    """Returns pointid of point on polydata closest to refpoint."""
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(polydata)
    locator.BuildLocator()
    return locator.FindClosestPoint(refpoint)


def gradient(polydata, inputarray, outputarray, iscelldata=False):
    """Compute the gradient of a scalar or vector field."""
    gradients = vtk.vtkGradientFilter()
    gradients.SetInput(polydata)
    fieldassociation = 1 if iscelldata == True else 0
    gradients.SetInputScalars(fieldassociation, inputarray)
    gradients.SetResultArrayName(outputarray)
    gradients.Update()
    return gradients.GetOutput()


def initializearray(polydata, arrayname, isscalar=True, ispointdata=True):
    """Initialize a data array. Choose pointdata or celldata and scalar or
    vector. Array values are initialized with 0.0 or [0.0, 0.0, 0.0]."""
    if ispointdata:  # add array to pointdata
        numberofpoints = polydata.GetNumberOfPoints()
        array = vtk.vtkDoubleArray()
        array.SetName(arrayname)
        if isscalar:  # array holds scalars
            array.SetNumberOfValues(numberofpoints)
            array.FillComponent(0, 0.0)
            polydata.GetPointData().AddArray(array)
        else:  # array holds vectors
            array.SetNumberOfComponents(3)
            array.SetNumberOfTuples(numberofpoints)
            for j in range(3):
                array.FillComponent(j, 0.0)
            polydata.GetPointData().AddArray(array)
    else:  # add array to celldata
        numberofcells = polydata.GetNumberOfCells()
        array = vtk.vtkDoubleArray()
        array.SetName(arrayname)
        if isscalar:  # array holds scalars
            array.SetNumberOfValues(numberofcells)
            array.FillComponent(0, 0.0)
            polydata.GetCellData().AddArray(array)
        else:  # array holds vectors
            array.SetNumberOfComponents(3)
            array.SetNumberOfTuples(numberofcells)
            for j in range(3):
                array.FillComponent(j, 0.0)
            polydata.GetCellData().AddArray(array)
    return array


def multiplyarray(polydata, inputarray, outputarray, value):
    """Multiply inputarray of polydata by specified value and store in
    outputarray."""
    calc = vtk.vtkArrayCalculator()
    calc.SetInput(polydata)
    calc.AddScalarArrayName(inputarray, 0)
    calc.SetFunction('{0} * {1}'.format(inputarray, value))
    calc.SetResultArrayName(outputarray)
    calc.Update()
    return calc.GetOutput()


def normalextrusion(polydata, extrusionlength=0.0, normalarrayname='Normals'):
    """Extrude points a specified length along the normals"""
    numberofpoints = polydata.GetNumberOfPoints()
    normalarray = polydata.GetPointData().GetArray(normalarrayname)
    newpoints = vtk.vtkPoints()
    newpoints.SetNumberOfPoints(numberofpoints)
    for i in range(numberofpoints):
        normal = [normalarray.GetComponent(i, 0),
                  normalarray.GetComponent(i, 1),
                  normalarray.GetComponent(i, 2)]
        point = polydata.GetPoint(i)
        newpoint = [point[0] + extrusionlength * normal[0],
                    point[1] + extrusionlength * normal[1],
                    point[2] + extrusionlength * normal[2]]
        if i == 0:
            print point, normal, newpoint
        newpoints.SetPoint(i, newpoint)
    polydata.SetPoints(newpoints)
    return polydata


def planeclip(polydata, point, normal, insideout=True):
    """Clip polydata with a plane defined by point and normal. Change clipping
    direction with 'insideout' argument."""
    clipplane = vtk.vtkPlane()
    clipplane.SetOrigin(point)
    clipplane.SetNormal(normal)
    clipper = vtk.vtkClipPolyData()
    clipper.SetInput(polydata)
    clipper.SetClipFunction(clipplane)
    if insideout:
        clipper.InsideOutOn()
    clipper.Update()
    return clipper.GetOutput()


def pointdatatocelldata(polydata):
    """Convert pointdata to celldata."""
    converter = vtk.vtkPointDataToCellData()
    converter.SetInput(polydata)
    converter.Update()
    return converter.GetOutput()


def pointnormals(surface):
    """Add point normals."""
    normals = vtk.vtkPolyDataNormals()
    normals.SetInput(surface)
    normals.ComputePointNormalsOn()
    normals.ComputeCellNormalsOff()
    normals.SplittingOff()
    normals.Update()
    return normals.GetOutput()


def probe(source, probe):
    """Compute point attributes (e.g. scalars, vectors, etc.) at all points of
    the probe object by interpolating the source data."""
    prober = vtk.vtkProbeFilter()
    prober.SetInput(probe)
    prober.SetSource(source)
    prober.Update()
    return prober.GetOutput()


def renderer(objectdicts, label=''):
    """Render several polydata in window.

    Args:
        objectdicts:
            List of Python dictionaries; one dictionary per polydata.
            The dictionaries can have three keys:
            1. polydata (required key; if not provided, dict will be skipped)
            2. color (optional; provide color in RGB (range: 0 - 1);
                default is [0, 0, 0])
            3. opacity (optional; 0 is fully transparent, 1 is fully opaque;
                default is 1)
        label: String to annotate the render window.

    Returns:
        vtkRenderWindow displaying polydata.

    Example usage:
        renderer([dict(polydata=surface, opacity=0.3),
                  dict(polydata=centerline, color=[1, 1, 1])],
                 label='All together!')

    TODO:
        * Close window with key press.
        * Place window in foreground.
        * When calling function in for-loop, window doesn't seem properly closed
          at end of iteration.

    """
    # Initialize render window
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)
    iren.SetRenderWindow(renWin)

    # Add each actor
    for objectdict in objectdicts:

        # If dict item exists use specified value, otherwise use default value
        polydata = objectdict.get('polydata')  # default is None
        color = objectdict.get('color', [1, 1, 1])
        opacity = objectdict.get('opacity', 1)

        if not polydata:  # only add actor if polydata is specified
            continue

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(polydata)
        mapper.ScalarVisibilityOff()
        actor = vtk.vtkActor()
        actor.GetProperty().SetColor(color)
        actor.GetProperty().SetOpacity(opacity)
        actor.SetMapper(mapper)
        ren.AddActor(actor)

    # Add text actor
    if label:
        txt = vtk.vtkTextActor()
        txt.SetInput(label)
        txtprop=txt.GetTextProperty()
        txtprop.SetFontFamilyToArial()
        txtprop.SetFontSize(18)
        txtprop.SetColor(0, 0, 0)
        txt.SetDisplayPosition(20, 30)
        ren.AddActor(txt)

    # Renderer (window) specifications
    ren.SetBackground(82/255., 87/255., 110/255.)  # Paraview default
    renWin.SetSize(640, 480)
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1)

    # Start rendering
    iren.Initialize()
    renWin.Render()
    iren.Start()


def scalarclip(polydata, arrayname, scalarvalue, insideout=True,
               ispointdata=True):
    """Clip vtkPolyData returning regions with value of array above a
    specified value."""
    clipper = vtk.vtkClipPolyData()
    clipper.SetInput(polydata)
    if ispointdata:  # array is pointdata
        clipper.SetInputArrayToProcess(0, 0, 0,
                                       vtk.vtkDataObject.
                                       FIELD_ASSOCIATION_POINTS,
                                       arrayname)
    else:  # array is celldata
        clipper.SetInputArrayToProcess(0, 0, 0,
                                       vtk.vtkDataObject.
                                       FIELD_ASSOCIATION_CELLS,
                                       arrayname)
    clipper.SetValue(scalarvalue)
    if insideout:  # inverse clipping direction
        clipper.InsideOutOn()
    clipper.Update()
    return clipper.GetOutput()


def slicedataset(dataset, point, normal):
    """Slice through a vtkDataSet object with a plane defined by point and
    normal."""
    cutplane = vtk.vtkPlane()
    cutplane.SetOrigin(point)
    cutplane.SetNormal(normal)
    cutter = vtk.vtkCutter()
    cutter.SetInput(dataset)
    cutter.SetCutFunction(cutplane)
    cutter.Update()
    return cutter.GetOutput()


def surfacearea(surface):
    """Compute surface area of input surface."""
    properties = vtk.vtkMassProperties()
    properties.SetInput(surface)
    properties.Update()
    return properties.GetSurfaceArea()


def surfacedistance(surface, referencesurface, distancearrayname='distance'):
    """Calculate unsigned distance to the reference surface at each point of the
    input surface."""
    distance = vtk.vtkDistancePolyDataFilter()
    distance.SetInput(0, surface)
    distance.SetInput(1, referencesurface)
    distance.SignedDistanceOff()
    distance.Update()
    osurface = distance.GetOutput()
    osurface.GetPointData().GetArray('Distance').SetName(distancearrayname)
    osurface.GetCellData().RemoveArray('Distance')
    return osurface


def surfacescaling(polydata, scalefactor=1.0):
    """Scale a vtkPolyData object."""
    transform = vtk.vtkTransform()
    transform.Scale(scalefactor, scalefactor, scalefactor)
    transformFilter = vtk.vtkTransformFilter()
    transformFilter.SetInput(polydata)
    transformFilter.SetTransform(transform)
    transformFilter.Update()
    return transformFilter.GetOutput()


def threshold(polydata, arrayname, valuerange=[0, 1], iscelldata=True):
    """Extract those cells from polydata whose pointdata/celldata values are
    within a specified range. For pointdata, cells are included if values of all
    cell points are within the range."""
    thresholder = vtk.vtkThreshold()
    thresholder.SetInput(polydata)
    thresholder.ThresholdBetween(valuerange[0], valuerange[1])
    if iscelldata:
        thresholder.SetInputArrayToProcess(0, 0, 0,
                                           vtk.vtkDataObject.
                                           FIELD_ASSOCIATION_CELLS, arrayname)
    else:
        thresholder.SetInputArrayToProcess(0, 0, 0,
                                           vtk.vtkDataObject.
                                           FIELD_ASSOCIATION_POINTS, arrayname)
    thresholder.Update()
    surfer = vtk.vtkDataSetSurfaceFilter()
    surfer.SetInput(thresholder.GetOutput())
    surfer.Update()
    return surfer.GetOutput()


def triangulate(surface):
    """Triangulate a surface mesh."""
    trianglefilter = vtk.vtkTriangleFilter()
    trianglefilter.SetInput(surface)
    trianglefilter.Update()
    return trianglefilter.GetOutput()


def vectormagnitude(surface, inputarray, outputarray, iscelldata=False):
    """Add data array with the magnitude of a vector field."""
    calc = vtk.vtkArrayCalculator()
    calc.SetInput(surface)
    if iscelldata:
        calc.SetAttributeModeToUseCellData()
    calc.AddVectorArrayName(inputarray, 0, 1, 2)
    calc.SetFunction('mag(%s)' % inputarray)
    calc.SetResultArrayName(outputarray)
    calc.Update()
    return calc.GetOutput()


def vectornormalization(surface, inputarray, outputarray, iscelldata=False):
    """Add data array with the normalization of a vector field."""
    calc = vtk.vtkArrayCalculator()
    calc.SetInput(surface)
    if iscelldata:
        calc.SetAttributeModeToUseCellData()
    calc.AddVectorArrayName(inputarray, 0, 1, 2)
    calc.SetFunction('norm(%s)' % inputarray)
    calc.SetResultArrayName(outputarray)
    calc.Update()
    return calc.GetOutput()


def volume(surface):
    """Compute volume of a closed triangulated surface mesh."""
    properties = vtk.vtkMassProperties()
    properties.SetInput(surface)
    properties.Update()
    return properties.GetVolume()


def vtkerror(path):
    """Write VTK error window text to file at specified path."""
    err = vtk.vtkFileOutputWindow()
    err.SetFileName(path)
    err.SetInstance(err)
