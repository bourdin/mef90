import ensight
import pymef90
import ensmef90
import os

e2cfile  = sys.argv[1]
#e2cfile  = '/Users/blaise/Desktop/DryingCylinder-R=2-L=2-h=0.1-GC=1.0e-2-DeltaT=10-24.old/DryingCylinder-R2-L2-h0.1.e2c'
genfile  = e2cfile.replace('e2c', 'gen')
enerfile = e2cfile.replace('e2c', 'ener')

laststep = pymef90.energetlaststep(enerfile)

### Computation specific parameters
NS_PartID=[2,4]
SeqMesh_PartID=[3]
DistMesh_PartID=[1]


### Read File
ensight.part.select_default()
ensight.part.modify_begin()
ensight.part.elt_representation("3D_border_2D_full")
ensight.part.modify_end()
ensight.data.binary_files_are("big_endian")
ensight.data.format("MultiExodusII")
ensight.data.reader_option("'Use distribution factors' ON")
ensight.data.reader_option("'Ignore side sets' OFF")
ensight.data.reader_option("'Ignore node sets' OFF")
ensight.data.reader_option("'Use node and element maps' ON")
ensight.data.reader_option("'Verbose mode' OFF")
ensight.data.reader_option("'Use higher-order elements' ON")
ensight.data.reader_option("'Ignore constant variables' OFF")
ensight.data.reader_option("'NaN filter input data' ON")
ensight.data.reader_option("'Clip overlapping timesteps' OFF")
ensight.data.reader_option("'Autodetect spatial decomp' OFF")
ensight.data.reader_option("'Use undef value for missing vars' ON")
ensight.data.reader_option("'Ignore elment attribute vars' OFF")
ensight.data.reader_option("'Use detected DTA XML file' ON")
ensight.data.reader_option("'Autogenerate DTA XML file' OFF")
ensight.data.reader_option("'Use full object names' OFF")
ensight.data.reader_option("'Epsilon' 1.000000e+00")
ensight.data.reader_option("'Scale factor' 1.000000e+00")
ensight.data.shift_time(1.0,0.0,0.0)
ensight.data.replace(e2cfile)

ensight.case.create_viewport("OFF")
ensight.case.apply_context("OFF")
ensight.case.reflect_model_in("'none'")
ensight.case.add("Geometry")
ensight.case.select("Geometry")
ensight.data.binary_files_are("big_endian")
ensight.data.format("MultiExodusIIng")
ensight.data.shift_time(1.0,0.0,0.0)
ensight.data.geometry(genfile)
ensight.data.read_all()

### Jump to last computed step
ensight.solution_time.current_step(laststep-1)
ensight.solution_time.update_to_current()


ensmef90.Init()
ensmef90.FractureActivate()

### hide node sets
ensight.part.select_begin( NS_PartID )
ensight.part.modify_begin()
ensight.part.visible("OFF")
ensight.part.modify_end()


### rotate everything
ensight.view_transf.rotate(1.702325e+01, -1.782179e+01, 0.000000e+00)
ensight.view_transf.fit(0)


### Make the Sequential mesh trasparent
ensight.part.select_begin(SeqMesh_PartID)
ensight.part.select_begin(3,4)
ensight.part.modify_begin()
ensight.part.colorby_rgb(0.6,0.6,0.6)
ensight.part.opaqueness(0.5)
ensight.part.modify_end()

### Extract Crack
ensmef90.CrackExtract(DistMesh_PartID, .05)
ensmef90.BrittlePart(5)

ensight.legend.select_palette_begin("_Fracture")
ensight.legend.visible("OFF")
ensight.function.palette("_Fracture")


### export PNG
ensmef90.ExportPNG(e2cfile)
