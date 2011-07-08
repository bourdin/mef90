import paraview.simple as pv

def ExodusIIReader(prefix,numproc=1,pattern='%s-%04i.gen'):
    filename = [pattern%(prefix,k) for k in range(numproc)]
    reader = pv.ExodusIIReader(FileName = filename)
    reader.FilePrefix=prefix
    reader.FilePattern = pattern
    reader.FileRange = 0,numproc-1
    reader.XMLFileName = ''
    reader.PointVariables = ['Displacement ','Fracture','Temperature']
    reader.ApplyDisplacements = 0
    return reader

def GetCrack(cutoff=.05,source=None):
    if not source:
        source = pv.GetActiveSource()

    Crack = pv.Contour( PointMergeMethod="Uniform Binning" )

    Crack.PointMergeMethod = "Uniform Binning"
    Crack.ContourBy = ['POINTS', 'Fracture']
    Crack.Isosurfaces = [cutoff]

    Crack.ComputeNormals = 0
    return Crack
    
#RenderView1 = GetRenderView()
#DataRepresentation1 = GetDisplayProperties(DryingBrickX1Y1Z1h_100)
#DataRepresentation2 = Show()
#DataRepresentation2.EdgeColor = [0.0, 0.0, 0.50000762951094835]