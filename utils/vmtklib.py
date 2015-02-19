"""Module for functions in the Vascular Modeling Toolkit (VMTK).

For more details on VMTK, go to www.vmtk.org. For more details specifically on
VMTK scripts, go to www.vmtk.org/vmtkscripts.

To a large extent, function names were chosen to match the names of
corresponding vmtk scripts. Parts of the documentation were taken from
www.vmtk.org.

"""

from vmtk import pypes, vmtkscripts


def vmtkbifurcationreferencesystems(centerlines):
    """Compute reference system for each bifurcation of a vessel tree.

    Args:
        centerlines: Centerlines split into branches.

    Returns:
        Reference system for each bifurcation.
        Pointdata (selection):
            Normal: Normal of the bifurcation plane.
            UpNormal: Normal pointing toward the bifurcation apex.

    """
    bifrefsystem = vmtkscripts.vmtkBifurcationReferenceSystems()
    bifrefsystem.Centerlines = centerlines
    bifrefsystem.RadiusArrayName = 'MaximumInscribedSphereRadius'
    bifrefsystem.BlankingArrayName = 'Blanking'
    bifrefsystem.GroupIdsArrayName = 'GroupIds'
    bifrefsystem.Execute()
    return bifrefsystem.ReferenceSystems


def vmtkbifurcationsections(surface, centerlines, distance=1):
    """Compute branch sections located a fixed number of maximally inscribed
    spheres away from each bifurcation.

    Args:
        surface: Surface split into branches.
        centerlines: Centerlines split into branches.
        distance: Distance from bifurcation in number of maximally inscribed
            spheres, where each sphere touches the center of the previous one.

    Returns:
        Polydata with one cross section per branch of each bifurcation.
        Celldata (selection):
            BifurcationSectionArea: Section area.
            BifurcationSectionMinSize: Minimum diameter of section.
            BifurcationSectionMaxSize: Maximum diameter of section.

    """
    bifsections = vmtkscripts.vmtkBifurcationSections()
    bifsections.Surface = surface
    bifsections.Centerlines = centerlines
    bifsections.NumberOfDistanceSpheres = distance
    bifsections.RadiusArrayName = 'MaximumInscribedSphereRadius'
    bifsections.GroupIdsArrayName = 'GroupIds'
    bifsections.CenterlineIdsArrayName = 'CenterlineIds'
    bifsections.TractIdsArrayName = 'TractIds'
    bifsections.BlankingArrayName = 'Blanking'
    bifsections.Execute()
    return bifsections.BifurcationSections


def vmtkbifurcationvectors(centerlines, referencesystems):
    """Compute bifurcation vectors.

    Args:
        centerlines: Centerlines split into branches.
        referencesystems: Bifurcation reference systems.

    Returns:
        Vectors in the direction of the parent and daughter branches of each
        bifurcation.
        Pointdata (selection):
            BifurcationVectors: Bifurcation vectors.
            InPlaneBifurcationVectors: BifurcationVectors projected onto the
                bifurcation plane.
            InPlaneBifurcationVectorAngles: Angles (in radians) between
                InPlaneBifurcationVectors and the UpNormal.

    """
    bifvectors = vmtkscripts.vmtkBifurcationVectors()
    bifvectors.Centerlines = centerlines
    bifvectors.ReferenceSystems = referencesystems
    bifvectors.RadiusArrayName = 'MaximumInscribedSphereRadius'
    bifvectors.GroupIdsArrayName = 'GroupIds'
    bifvectors.CenterlineIdsArrayName = 'CenterlineIds'
    bifvectors.TractIdsArrayName = 'TractIds'
    bifvectors.BlankingArrayName = 'Blanking'
    bifvectors.ReferenceSystemsNormalArrayName = 'Normal'
    bifvectors.ReferenceSystemsUpNormalArrayName = 'UpNormal'
    bifvectors.Execute()
    return bifvectors.BifurcationVectors


def vmtkbranchclipper(surface, centerlines, groupidlist=[], insideout=0):
    """Split surface into branches.

    Args:
        surface: Surface mesh of vascular geometry.
        centerlines: Centerlines split into branches.
        groupidlist: List of branches (specified by their GroupIds) to extract.
        insideout (bool): Either keep or remove the branches in groupidlist.

    Returns:
        Surface split into branches or, if specified, only the selected
        branches.

    """
    clipper = vmtkscripts.vmtkBranchClipper()
    clipper.Surface = surface
    clipper.Centerlines = centerlines
    clipper.InsideOut = insideout
    clipper.GroupIds = groupidlist
    clipper.GroupIdsArrayName = 'GroupIds'
    clipper.RadiusArrayName = 'MaximumInscribedSphereRadius'
    clipper.BlankingArrayName = 'Blanking'
    clipper.Execute()
    return clipper.Surface


def vmtkbranchextractor(centerlines):
    """Split centerlines into branches.

    Args:
        centerlines: Centerlines.

    Returns:
        Centerlines split into branches.
        Celldata (selection):
            CenterlineId: Cellid of centerline from which the tract was split.
            TractId: Id identifying each tract along one centerline.
            GroupId: Id of the group to which the tract belongs.
            Blanking: Boolean indicating whether tract belongs to bifurcation.

    """
    extractor = vmtkscripts.vmtkBranchExtractor()
    extractor.Centerlines = centerlines
    extractor.RadiusArrayName = 'MaximumInscribedSphereRadius'
    extractor.Execute()
    return extractor.Centerlines


def vmtkbranchgeometry(centerlines, smoothing=0, iterations=100):
    """Compute geometric variables for each branch of a vessel tree.

    Args:
        centerlines: Centerlines split into branches.
        smoothing (bool): Laplacian smooth branches before computing geometric
            variables.
        iterations: Number of smoothing iterations.

    Returns:
        Description of geometry for each branch. Vertices at (0, 0, 0), one per
        branch, are used as placeholders to assign pointdata with values
        for the geometric variables.
        Pointdata (selection):
            Length: Branch length.
            Curvature: Average curvature of branch.
            Torsion: Average torsion of branch.
            Tortuosity: Tortuosity of branch.

    Note:
        Smoothing doesn't seem to work, it does work with
        vmtkcenterlinegeometry.

    """
    branchgeometry = vmtkscripts.vmtkBranchGeometry()
    branchgeometry.Centerlines = centerlines
    branchgeometry.GroupIdsArrayName = 'GroupIds'
    branchgeometry.RadiusArrayName = 'MaximumInscribedSphereRadius'
    branchgeometry.BlankingArrayName = 'Blanking'
    branchgeometry.LineSmoothing = smoothing
    branchgeometry.NumberOfSmoothingIterations = iterations
    branchgeometry.Execute()
    return branchgeometry.GeometryData


def vmtkbranchmapping(surface, centerlines, referencesystems):
    """Map and stretch the longitudinal metrics obtained with
    vmtkbranchmetrics.

    Args:
        surface: Surface split into branches.
        centerlines: Centerlines split into branches.
        referencesystems: Bifurcation reference systems.

    Returns:
        Surface with the longitudinal metric mapped and stretched to correctly
        account for the presence of insertion regions at bifurcations.
        Pointdata (selection):
            StretchedMapping: Corrected longitudinal metric.

    """
    mapper = vmtkscripts.vmtkBranchMapping()
    mapper.Surface = surface
    mapper.Centerlines = centerlines
    mapper.ReferenceSystems = referencesystems
    mapper.AbscissasArrayName = 'Abscissas'
    mapper.NormalsArrayName = 'ParallelTransportNormals'
    mapper.GroupIdsArrayName = 'GroupIds'
    mapper.CenterlineIdsArrayName = 'CenterlineIds'
    mapper.TractIdsArrayName = 'TractIds'
    mapper.ReferenceSystemsNormalArrayName = 'Normal'
    mapper.RadiusArrayName = 'MaximumInscribedSphereRadius'
    mapper.BlankingArrayName = 'Blanking'
    mapper.AngularMetricArrayName = 'AngularMetric'
    mapper.AbscissaMetricArrayName = 'AbscissaMetric'
    mapper.Execute()
    return mapper.Surface


def vmtkbranchmetrics(surface, centerlines):
    """Compute longitudinal and circumferential metrics for each branch.

    Args:
        surface: Surface split into branches.
        centerlines: Centerlines split into branches.

    Returns:
        Surface with longitudinal and circumferential metrics for each
        branch.
        Pointdata (selection):
            AbscissaMetric: Curvilinear abscissa of centerlines projected onto
                the surface.
            AngularMetric: Periodic circumferential coordinates of surface mesh
                points around the centerlines, spanning the interval (-pi, pi).
                The zero angle is derived from the ParallelTransportNormals
                pointdata of the centerlines, which can be obtained using
                vmtkcenterlineattributes.

    """
    branchmetrics = vmtkscripts.vmtkBranchMetrics()
    branchmetrics.Surface = surface
    branchmetrics.Centerlines = centerlines
    branchmetrics.AbscissasArrayName = 'Abscissas'
    branchmetrics.NormalsArrayName = 'ParallelTransportNormals'
    branchmetrics.RadiusArrayName = 'MaximumInscribedSphereRadius'
    branchmetrics.GroupIdsArrayName = 'GroupIds'
    branchmetrics.CenterlineIdsArrayName = 'CenterlineIds'
    branchmetrics.TractIdsArrayName = 'TractIds'
    branchmetrics.BlankingArrayName = 'Blanking'
    branchmetrics.Execute()
    return branchmetrics.Surface


def vmtkbranchpatching(surface, longitudinalpatchsize=1.0,
                       circularnumberofpatches=12):
    """Patch surface of each branch.

    Patching means to 'cut' a set of contiguous rectangular regions on the
    surface mesh that follow iso-contours in the StretchedMapping and
    AngularMetric arrays. All the quantities of interest (e.g. wall shear
    stress, oscillatory shear index) are averaged on each of these patches.

    Args:
        surface: Surface split into branches.
        longitudinalpatchsize: 'Length' of the patch along the longitudinal
            direction.
        circularnumberofpatches: Number of patches along the circumference.

    Returns:
        Surface composed of disconnected patches.
        Pointdata:
            All the original pointdata
        Celldata (selection):
            Average of original pointdata on patch to which cell belongs
            Slab: Patch coordinate (integer) in longitudinal direction.
            Sector: Patch coordinate (integer) in circumferential direction.
            PatchArea: Surface area of the patch

    """
    patcher = vmtkscripts.vmtkBranchPatching()
    patcher.Surface = surface
    patcher.LongitudinalPatchSize = longitudinalpatchsize
    patcher.CircularNumberOfPatches = circularnumberofpatches
    patcher.UseConnectivity = 1
    patcher.CircularPatching = 1
    patcher.GroupIdsArrayName = 'GroupIds'
    patcher.LongitudinalMappingArrayName = 'StretchedMapping'
    patcher.CircularMappingArrayName = 'AngularMetric'
    patcher.Execute()
    #return (patcher.Surface, patcher.PatchedData)
    return patcher.Surface


def vmtkcenterlineattributes(centerlines):
    """Compute centerline attributes.

    Args:
        centerlines: Centerlines.

    Returns:
        Centerlines with attributes.
        Pointdata:
            MaximumInscribedSphereRadius: If the point on the centerline is the
                center of a sphere, this is the radius of the largest possible
                sphere that does not intersect the surface.
            Abscissas: Position along the centerlines. By default, the abscissa
                is measured from the start of the centerlines.
            ParallelTransportNormals: 'Normal' of the centerlines (perpendicular
                to centerline direction).

    """
    clattributes = vmtkscripts.vmtkCenterlineAttributes()
    clattributes.Centerlines = centerlines
    clattributes.Execute()
    return clattributes.Centerlines


def vmtkcenterlinegeometry(centerlines, smoothing=0, iterations=100):
    """Compute the local geometry of centerlines.

    Args:
        centerlines: Centerlines.
        smoothing (bool): Laplacian smooth centerlines before computing
            geometric variables.
        iterations: Number of smoothing iterations.

    Returns:
        Centerlines with geometric variables defined at each point.
        Pointdata (selection):
            Curvature: Local curvature.
            Torsion: Local torsion.
        Celldata (selection):
            Tortuosity: Tortuosity of each centerline.
            Length: Length of each centerline.

    Note:
        Since the computation of the geometric variables depends on first,
        second and third derivatives of the line coordinates, and since such
        derivatives are approximated using a simple finite difference scheme
        along the line, it is very likely that such derivatives will be affected
        by noise that is not appreciable when looking at the line itself. For
        this reason, it might be necessary to run the Laplacian smoothing filter
        before computing the derivatives and the related quantities.

    """
    clgeometry = vmtkscripts.vmtkCenterlineGeometry()
    clgeometry.Centerlines = centerlines
    clgeometry.LineSmoothing = smoothing
    clgeometry.NumberOfSmoothingIterations = iterations
    clgeometry.Execute()
    return clgeometry.Centerlines


def vmtkcenterlinemerge(centerlines, length=.1):
    """Merge centerline tracts belonging to the same groups.

    Args:
        centerlines: Centerlines split into branches.
        length: Distance between centerline points after resampling.

    Returns:
        Centerlines with only one centerline branch per vessel tree branch. The
        centerline branches meet at the bifurcation origins.

    """
    merger = vmtkscripts.vmtkCenterlineMerge()
    merger.Centerlines = centerlines
    merger.Length = length
    merger.RadiusArrayName = 'MaximumInscribedSphereRadius'
    merger.GroupIdsArrayName = 'GroupIds'
    merger.CenterlineIdsArrayName = 'CenterlineIds'
    merger.BlankingArrayName = 'Blanking'
    merger.TractIdsArrayName = 'TractIds'
    merger.Execute()
    return merger.Centerlines


def vmtkcenterlinemodeller(centerlines, size=[64, 64, 64]):
    """Convert a centerline to an image containing the tube function.

    Args:
        centerlines: Centerlines.
        size: Image dimensions.

    Returns:
        Signed distance transform image, with the zero level set being (tapered)
        tubes running from one centerline point to the next with a radius at
        each end corresponding to the local MaximumInscribedSphereRadius.

    """
    modeller = vmtkscripts.vmtkCenterlineModeller()
    modeller.Centerlines = centerlines
    modeller.RadiusArrayName = 'MaximumInscribedSphereRadius'
    modeller.SampleDimensions = size
    modeller.Execute()
    return modeller.Image


def vmtkcenterlineoffsetattributes(centerlines, referencesystems,
    refgroupid=1):
    """Offset centerline attributes to a bifurcation reference system.

    Args:
        centerlines: Centerlines with attributes.
        referencesystems: Bifurcation reference systems.
        refgroupid: GroupId of bifurcation to which to offset the attributes.

    Returns:
        Centerlines with the attributes offset in such a way that the abscissa
        of the closest point to the bifurcation origin is zero and the
        centerline normal at that point coincides with the bifurcation reference
        system normal.

    """
    offsetter = vmtkscripts.vmtkCenterlineOffsetAttributes()
    offsetter.Centerlines = centerlines
    offsetter.ReferenceSystems = referencesystems
    offsetter.ReferenceGroupId = refgroupid
    offsetter.AbscissasArrayName = 'Abscissas'
    offsetter.NormalsArrayName = 'ParallelTransportNormals'
    offsetter.GroupIdsArrayName = 'GroupIds'
    offsetter.CenterlineIdsArrayName = 'CenterlineIds'
    offsetter.ReferenceSystemsNormalArrayName = 'Normal'
    offsetter.Execute()
    return offsetter.Centerlines


def vmtkcenterlineresampling(centerlines, length=.1):
    """Resample input centerlines with a spline filter.

    Args:
        centerlines: Centerlines.
        length: Space between centerline points after resampling.

    Returns:
        Centerlines with equal spacing between centerline points.

    """
    resampler = vmtkscripts.vmtkCenterlineResampling()
    resampler.Centerlines = centerlines
    resampler.Length = length
    resampler.Execute()
    return resampler.Centerlines


def vmtkcenterlines(surface, endpoints=0, interactive=False,
                    sourcepoints=[0, 0, 0], targetpoints=[0, 0, 0]):
    """Compute centerlines of a vascular geometry.

    Args:
        surface: Surface mesh of a vascular geometry.
        endpoints (bool): Include endpoints. By construction, centerlines do
            not reach the source/targetpoints. endpoints=1 bridges the
            start/end of the centerlines to the source/targetpoints.
        interactive (bool): Select source/targetpoints interactively. This
            pops up a VTK window. Follow instructions in terminal.
        sourcepoints: Give barycenters of sourcepoints as in '[0.0, 1.0, 2.0]'.
            In case of multiple sourcepoints, append the lists of each three
            coordinates as in '[0.0, 1.0, 2.0, 3.0, 4.0, 5.0]'.
        targetpoints: Give barycenters of targetpoints following same notation
            as with sourcepoints.

    Returns:
        Centerlines. Each cell of the centerlines polydata is a centerline from
        one of the sourcepoints to one of the targetpoints.

    """
    centerliner = vmtkscripts.vmtkCenterlines()
    centerliner.Surface = surface
    centerliner.AppendEndPoints = endpoints
    if interactive:
        centerliner.SeedSelectorName = 'pickpoint'
    else:
        centerliner.SeedSelectorName = 'pointlist'
        centerliner.SourcePoints = sourcepoints
        centerliner.TargetPoints = targetpoints
    centerliner.Execute()
    return centerliner.Centerlines


def vmtkcenterlinesections(surface, centerlines):
    """Compute sections located at each point of the centerlines.

    Args:
        surface: Surface mesh of vascular geometry.
        centerlines: Centerlines corresponding to surface.

    Returns:
        Polydata with one cross section per branch of each bifurcation.
        Celldata (selection):
            CenterlineSectionArea: Section area.
            CenterlineSectionMinSize: Minimum diameter of section.
            CenterlineSectionMaxSize: Maximum diameter of section.

    """
    sectioner = vmtkscripts.vmtkCenterlineSections()
    sectioner.Surface = surface
    sectioner.Centerlines = centerlines
    sectioner.Execute()
    return sectioner.CenterlineSections


def vmtkcenterlinesmoothing(centerlines, iterations=100):
    """Smooth centerlines with a Laplacian smoothing filter.

    Args:
        centerlines: Centerlines.
        iterations: Number of smoothing iterations.

    Returns:
        Smoothed centerlines.

    """
    smoother = vmtkscripts.vmtkCenterlineSmoothing()
    smoother.Centerlines = centerlines
    smoother.NumberOfSmoothingIterations = iterations
    smoother.Execute()
    return smoother.Centerlines


def vmtkdelaunayvoronoi(surface, removesubresolution=1):
    """Compute Voronoi diagram corresponding to the surface.

    Args:
        surface: Surface mesh.
        removesubresolution (bool): Remove Voronoi points with 'subresolution'
            radii.

    Returns:
        Voronoi diagram.

    Note:
        The Voronoi diagram represents the lumen volume as a set of overlapping
        spheres: the largest spheres correspond to the local minimum lumen
        diameter, the smallest spheres correspond to finer details on the
        surface.

    """
    voronoi = vmtkscripts.vmtkDelaunayVoronoi()
    voronoi.Surface = surface
    voronoi.RemoveSubresolutionTetrahedra = removesubresolution
    voronoi.Execute()
    return voronoi.VoronoiDiagram


def vmtkdistancetocenterlines(surface, centerlines):
    """Compute distance from the surface to the centerlines.

    Args:
        surface: Surface mesh of vascular geometry.
        centerlines: Centerlines corresponding to surface.

    Returns:
        Surface with DistanceToCenterlines as pointdata.

    """
    distance = vmtkscripts.vmtkDistanceToCenterlines()
    distance.Surface = surface
    distance.Centerlines = centerlines
    distance.RadiusArrayName = 'MaximumInscribedSphereRadius'
    distance.Execute()
    return distance.Surface


def vmtkflowextensions(surface, interactive=1, extensionlength=10,
                       transitionratio=.25):
    """Extrude inlets and outlets of a vascular geometry.

    Args:
        surface: Surface mesh of vascular geometry.
        interactive (bool): Choose inlets/outlets to be extruded.
        extensionlength: Length of extrusions.
        transitionratio: Rate of transition from model section to circular
        section.

    Returns:
        Surface with extruded inlets/outlets.

    """
    extender = vmtkscripts.vmtkFlowExtensions()
    extender.Surface = surface
    extender.ExtensionMode = 'boundarynormal'
    extender.TransitionRatio = transitionratio
    extender.Interactive = interactive
    extender.ExtensionLength = extensionlength
    extender.Execute()
    return extender.Surface


def vmtkicpregistration(surface, referencesurface):
    """Register a surface to a reference surface using the interative closest
    point algorithm.

    Args:
        surface = Surface mesh.
        referencesurface = Reference surface mesh.

    Returns:
        'surface' rigidly transformed to best match 'referencesurface'.
        Pointdata:
            Distance: Distance from 'surface' to 'referencesurface'.

    """
    registration = vmtkscripts.vmtkICPRegistration()
    registration.Surface = surface
    registration.ReferenceSurface = referencesurface
    registration.DistanceArrayName = 'Distance'
    registration.Execute()
    return registration.Surface


def vmtkimagereader(path):
    """Read an image and store it in a vtkImageData object.

    Args:
        path: Path to the image file.

    Returns:
        vtkImageData object.

    Note:
        Reads several image formats: vti, vtk, dcm, raw, mha, mhd, tif, png

    """
    reader = vmtkscripts.vmtkImageReader()
    reader.InputFileName = path
    reader.Execute()
    return reader.Image


def vmtkimagewriter(image, path):
    """Write a vtkImageData object to disk.

    Args:
        image: vtkImageData object.
        path: Path to the image file.

    Returns:
        n/a

    Note:
        Writes several image formats: vti, vtk, mha, mhd, tif, png, dat

    """
    writer = vmtkscripts.vmtkImageWriter()
    writer.Image = image
    writer.OutputFileName = path
    writer.Execute()


def vmtkmarchingcubes(image, level=0.0):
    """Generate an isosurface of given level from a 3D image.

    Args:
        image: vtkImageData object.
        level: Graylevel at which to generate the isosurface.

    Returns:
        Surface mesh of the isosurface.

    """
    marcher = vmtkscripts.vmtkMarchingCubes()
    marcher.Image = image
    marcher.Level = level
    marcher.Execute()
    return marcher.Surface


def vmtkmeshreader(path):
    """Read a mesh and store it in a vtkUnstructuredGrid object.

    Args:
        path: Path to the mesh file.

    Returns:
        vtkUnstructuredGrid object.

    Note:
        Reads several mesh formats: vtu, vtk, FDNEUT, xda, neu (ngneut),
        gneu (gambit), tec (tecplot), node (tetgen), ele (tetgen)

    """
    reader = vmtkscripts.vmtkMeshReader()
    reader.InputFileName = path
    reader.Execute()
    return reader.Mesh


def vmtkmeshtosurface(mesh, cleanoutput=1):
    """Convert a mesh to a surface by throwing out volume elements and (optionally) the relative points

    Args:
        mesh: Volumetric mesh.
        cleanoutput (bool): Remove unused points.

    Returns:
        vtkPolyData object.

    """
    extractor = vmtkscripts.vmtkMeshToSurface()
    extractor.Mesh = mesh
    extractor.CleanOutput = cleanoutput
    extractor.Execute()
    return extractor.Surface


def vmtkmeshvectorfromcomponents(mesh, vectorname='Velocity',
                                 componentsnames=['VelocityX', 'VelocityY',
                                 'VelocityZ'], removecomponents=False):
    """Create a vector array from a number of scalar arrays treated as vector
    components.

    Args:
        mesh: vtkUnstructuredGrid object.
        vectorname: Name to give to the vector array that will be created.
        componentsnames: List of names of the three scalar arrays that will be
            used as vector components.
        removecomponents (bool): Remove the scalar arrays after creating the
            vector array.

    Returns:
        vtkUnstructuredGrid object with vector array.

    """
    vectorer = vmtkscripts.vmtkMeshVectorFromComponents()
    vectorer.Mesh = mesh
    vectorer.VectorArrayName = vectorname
    vectorer.ComponentsArrayNames = componentsnames
    vectorer.RemoveComponentArrays = removecomponents
    vectorer.Execute()
    return vectorer.Mesh


def vmtkmeshwriter(mesh, path):
    """Write a vtkUnstructuredGrid object to disk.

    Args:
        image: vtkUnstructuredGrid object.
        path: Path to the mesh file.

    Returns:
        n/a

    Note:
        Writes several mesh formats: vtu, vtk, xda, FDNEUT, lifev, xml
        (dolfin), msh (fluent), tec (tecplot), node (tetgen), ele (tetgen), dat

    """
    writer = vmtkscripts.vmtkMeshWriter()
    writer.Mesh = mesh
    writer.OutputFileName = path
    writer.Execute()


def vmtknetworkextraction(surface):
    """Extract a network of approximated centerlines from a surface.

    Args:
        surface: Surface mesh of vascular network.

    Returns:
        Network of centerlines.

    Note:
        The surface must have at least one opening.

    """
    extractor = vmtkscripts.vmtkNetworkExtraction()
    extractor.Surface = surface
    extractor.Execute()
    return extractor.Network


def vmtkpolyballmodeller(voronoi, size=[64, 64, 64]):
    """Converts a polyball to an image containing the tube function.

    Args:
        voronoi: Voronoi diagram.
        size: Image dimensions.

    Returns:
        Signed distance transform image, with the zero level set being the lumen
        surface represented by the Voronoi diagram.

    Note:
        The Voronoi diagram represents the lumen volume as a set of overlapping
        spheres (also named 'polyball'): the largest spheres correspond to the
        local minimum lumen diameter, the smallest spheres correspond to finer
        details on the surface.

    """
    modeller = vmtkscripts.vmtkPolyBallModeller()
    modeller.Surface = voronoi
    modeller.RadiusArrayName = 'MaximumInscribedSphereRadius'
    modeller.SampleDimensions = size
    modeller.Execute()
    return modeller.Image


def vmtkpointsplitextractor(centerlines, splitpoint, gap=1.0):
    """Split centerlines at specified location.

    Args:
        centerlines: Centerlines.
        splitpoint: Location where to split the centerlines.
        gap: Length of 'Blanking=1' part of the centerlines.

    Returns
        Centerlines split at splitpoint, with the center of the gap at the
        splitpoint. The output is similar to the output of
        vmtkbranchextractor.
        Celldata (selection):
            CenterlineId: Cellid of centerline from which the tract was split.
            TractId: Id identifying each tract along one centerline.
            GroupId: Id of the group to which the tract belongs.
            Blanking: Boolean indicating whether tract belongs to bifurcation.

    """
    extractor = vmtkscripts.vmtkPointSplitExtractor()
    extractor.Centerlines = centerlines
    extractor.RadiusArrayName = 'MaximumInscribedSphereRadius'
    extractor.GroupIdsArrayName = 'GroupIds'
    extractor.SplitPoint = splitpoint
    extractor.Execute()
    return extractor.Centerlines


def vmtksurfacecapper(surface):
    """Caps the holes of a surface.

    Args:
        surface: Surface mesh of vascular geometry with holes at inlets and
            outlets.

    Returns:
        Surface mesh with capped holes. Each cap has an ID assigned for easy
        specification of boundary conditions.
        Celldata:
            CellEntityIds: ID assigned to caps.

    """
    capper = vmtkscripts.vmtkSurfaceCapper()
    capper.Surface = surface
    capper.Method = 'centerpoint'
    capper.Interactive = 0
    capper.Execute()
    return capper.Surface


def vmtksurfacecenterlineprojection(surface, centerlines):
    """Project pointdata from centerlines to a surface.

    Args:
        surface: Surface mesh of vascular geometry.
        centerlines: Centerlines corresponding to surface.

    Returns:
        Surface mesh with centerlines pointdata projected onto it.

    """
    projection = vmtkscripts.vmtkSurfaceCenterlineProjection()
    projection.Surface = surface
    projection.Centerlines = centerlines
    projection.RadiusArrayName = 'MaximumInscribedSphereRadius'
    projection.Execute()
    return projection.Surface


def vmtksurfacedecimation(surface, target=.9):
    """Reduce the number of triangles in a surface.

    Args:
        surface: Surface mesh.
        target: Desired number of triangles relative to the input number of
            triangles.

    Returns:
        Surface mesh with fewer triangles.

    """
    decimator = vmtkscripts.vmtkSurfaceDecimation()
    decimator.Surface = surface
    decimator.TargetReduction = target
    decimator.Execute()
    return decimator.Surface


def vmtksurfacedistance(surface, referencesurface,
                        distancearrayname='Distance'):
    """Compute the pointwise minimum distance from a surface to a reference
    surface.

    Args:
        surface: Surface mesh.
        referencesurface: Reference surface mesh.

    Returns:
        'surface' with 'Distance' pointdata.

    """
    distancer = vmtkscripts.vmtkSurfaceDistance()
    distancer.Surface = surface
    distancer.ReferenceSurface = referencesurface
    distancer.DistanceArrayName = distancearrayname
    distancer.DistanceVectorsArrayName = 'Distance'
    distancer.SignedDistanceArrayName = ''
    distancer.Execute()
    return distancer.Surface


def vmtksurfacenormals(surface):
    """Compute normals to a surface.

    Args:
        surface: Surface mesh.

    Returns:
        Surface mesh with 'Normals' vector pointdata.

    """
    normaller = vmtkscripts.vmtkSurfaceNormals()
    normaller.Surface = surface
    normaller.Execute()
    return normaller.Surface


def vmtksurfaceprojection(surface, referencesurface):
    """Interpolates the pointdata of a reference surface onto a surface based
    on minimum distance criterion.

    Args:
        surface: Surface mesh.
        referencesurface: Reference surface mesh.

    Returns:
        'surface' with projected pointdata from 'referencesurface'.

    """
    projector = vmtkscripts.vmtkSurfaceProjection()
    projector.Surface = surface
    projector.ReferenceSurface = referencesurface
    projector.Execute()
    return projector.Surface


def vmtksurfacereader(path):
    """Read a polydata (surface or centerline) and store it in a vtkPolyData
    object.

    Args:
        path: Path to the polydata file.

    Returns:
        vtkPolyData object.

    Note:
        Reads several polydata formats: vtp, vtk, stl, ply, tec (tecplot),
        dat (tecplot)

    """

    reader = vmtkscripts.vmtkSurfaceReader()
    reader.InputFileName = path
    reader.Execute()
    return reader.Surface


def vmtksurfaceremeshing(surface, edgelength=1.0, iterations=10):
    """Remesh a surface using high quality triangles.

    Args:
        surface: Surface mesh.
        edgelength: Target length of triangle edges.
        iterations: Number of iterations to optimize the mesh quality.

    Returns:
        Remeshed surface.

    """
    remesher = vmtkscripts.vmtkSurfaceRemeshing()
    remesher.Surface = surface
    remesher.NumberOfIterations = iterations
    remesher.ElementSizeMode = 'edgelength'
    remesher.TargetEdgeLength = edgelength
    remesher.Execute()
    return remesher.Surface


def vmtksurfacescaling(polydata, scalefactor):
    """Scale a polydata by an isotropic factor.

    Args:
        polydata: Surface mesh or centerlines.
        scalefactor: Scaling factor.

    Returns:
        Scaled polydata.

    """
    scaler = vmtkscripts.vmtkSurfaceScaling()
    scaler.Surface = polydata
    scaler.ScaleFactor = scalefactor
    scaler.Execute()
    return scaler.Surface


def vmtksurfacesmoothing(surface, iterations=100, method='taubin'):
    """Smooth a surface.

    Args:
        surface: Surface mesh.
        iterations: Number of smoothing iterations.
        method ('taubin', 'laplace'): Taubin's volume-preserving or a Laplacian
            smoothing filter.

    Returns:
        Smoothed surface.

    """
    smoother = vmtkscripts.vmtkSurfaceSmoothing()
    smoother.Surface = surface
    smoother.Method = method
    smoother.NumberOfIterations = iterations
    smoother.PassBand = 0.1
    smoother.Execute()
    return smoother.Surface


def vmtksurfacesubdivision(surface, numberofsubdivisions=1, method='linear'):
    """Subdivide a triangulated surface.

    Args:
        surface: Surface mesh.
        numberofsubdivisions: Number of subdivisions.
        method ('linear', 'butterfly', 'loop'): Subdivision method.

    Returns:
        Subdivided surface mesh.

    """
    divider = vmtkscripts.vmtkSurfaceSubdivision()
    divider.Surface = surface
    divider.NumberOfSubdivisions = numberofsubdivisions
    divider.Method = method
    divider.Execute()
    return divider.Surface


def vmtksurfacewriter(polydata, path):
    """Write a vtkPolyData object (e.g. surface and centerlines) to disk.

    Args:
        polydata: vtkPolyData object.
        path: Path to the polydata file.

    Returns:
        n/a

    Note:
        Writes several polydata formats: vtp, vtk, stl (use only for
        triangulated surface meshes), ply, tec (tecplot), dat

    """
    writer = vmtkscripts.vmtkSurfaceWriter()
    writer.Surface = polydata
    writer.OutputFileName = path
    writer.Execute()
