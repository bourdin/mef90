def init()
  ensight.data.binary_files_are("big_endian")
  ensight.command.part_selection_by("number")
  ensight.data.binary_files_are("big_endian")
  ensight.data.format("MultiExodusII")
  ensight.data.reader_option("'Ignore node sets' OFF")
  ensight.data.reader_option("'Use node and element maps' ON")
  ensight.data.reader_option("'Verbose mode' OFF")
  ensight.data.reader_option("'Ignore constant variables' OFF")
  ensight.data.reader_option("'NaN filter input data' ON")
  ensight.data.reader_option("'Clip overlapping timesteps' OFF")
  ensight.data.reader_option("'Autodetect spatial decomp' ON")
  ensight.data.reader_option("'Use undef value for missing vars' ON")
  ensight.data.reader_option("'Ignore elment attribute vars' OFF")
  ensight.data.reader_option("'Use detected DTA XML file' ON")
  ensight.data.reader_option("'Epsilon' 1.000000e+00")
  ensight.data.reader_option("'Scale factor' 1.000000e+00")
  ensight.data.shift_time(1.000000,0.000000,0.000000)
  
  ###
  ### Black Background
  ###
  ensight.viewport.select_begin(0)
  ensight.view_transf.function("global")
  ensight.viewport.background_type("constant")
  
  ###
  ### remove axis
  ###
  ensight.annotation.axis_model("OFF")

def FractureActivate():
  ensight.variables.activate("_Fracture")
  
  ###
  ### Set colormap for _Fracture
  ###
  ensight.legend.select_palette_begin("_Fracture")
  ensight.legend.visible("ON")
  ensight.function.palette("_Fracture")
  ensight.function.modify_begin()
  ensight.function.edit_level(1)
  ensight.function.rgb(1.0000e+00,0.0000e+00,0.0000e+00)
  ensight.function.edit_level(2)
  ensight.function.rgb(1.0000e+00,1.0000e+00,0.0000e+00)
  ensight.function.edit_level(3)
  ensight.function.rgb(0.0000e+00,1.0000e+00,0.0000e+00)
  ensight.function.edit_level(4)
  ensight.function.rgb(0.0000e+00,1.0000e+00,1.0000e+00)
  ensight.function.edit_level(5)
  ensight.function.rgb(0.0000e+00,0.0000e+00,1.0000e+00)
  ensight.function.range(0.000000,1.000000)
  ensight.function.modify_end()


def BrittlePart(partlist):
  ensight.part.select_begin(partlist)
  ensight.part.modify_begin()
  ensight.part.colorby_palette("_Fracture")
  ensight.part.modify_end()
  ensight.part.select_end()

def ElasticPart(partlist):
  ensight.part.select_begin(partlist)
  ensight.part.modify_begin()
  ensight.part.colorby_palette("_Fracture")
  ensight.part.opaqueness(4.000000e-01)
  ensight.part.modify_end()
  ensight.part.select_end()

def DispPart(partlist, factor):
  ensight.part.select_begin(partlist)
  ensight.part.select_end()
  ensight.part.modify_begin()
  ensight.part.displace_by("_Displacement_VEC")
  ensight.part.displace_factor(factor)
  ensight.part.modify_end()
  ensight.part.select_end()

#def ForcePart(partlist):

def CrackExtract(partlist, threshold):
  ensight.part.select_begin(partlist)
  ensight.part.modify_begin()
  ensight.part.visible("OFF")
  ensight.part.modify_end()
  ensight.part.modify_begin()
  ensight.isos.variable("_Fracture")
  ensight.isos.type("isovolume")
  ensight.isos.constraint("low")
  ensight.isos.min(threshold)
  ensight.isos.max(1.000000e+00)
  ensight.isos.create()
  ensight.part.modify_end()
  ensight.part.select_end()

def CrackRemove(partlist, threshold):
  ensight.part.select_begin(partlist)
  ensight.part.modify_begin()
  ensight.part.visible("OFF")
  ensight.part.modify_end()
  ensight.part.modify_begin()
  ensight.isos.variable("_Fracture")
  ensight.isos.type("isovolume")
  ensight.isos.constraint("high")
  ensight.isos.main(0.000000e+00)
  ensight.isos.max(threshold)
  ensight.isos.create()
  ensight.part.modify_end()
  ensight.part.select_end()

#def TempIsolines(partlist):

def ExportPNG(filename):
  ensight.file.image_file(filename)
  ensight.file.image_format("png")
  ensight.file.image_format_options("Compression Default")
  ensight.anim_recorders.render_offscreen("ON")
  ensight.file.image_numpasses(4)
  ensight.file.image_stereo("current")
  ensight.file.image_screen_tiling(1,1)
  ensight.file.image_window_size("HD1080p")
  ensight.file.save_image()

def ExportMOV(filename):
  ensight.file.animation_file(filename)
  ensight.file.image_format("mov")
  ensight.anim_recorders.render_offscreen("ON")
  ensight.file.image_numpasses(4)
  ensight.file.image_stereo("current")
  ensight.file.image_screen_tiling(1,1)
  ensight.file.image_file(filename)
  ensight.anim_recorders.render_offscreen("ON")
  ensight.file.image_numpasses(4)
  ensight.file.image_stereo("current")
  ensight.file.image_screen_tiling(1,1)
  ensight.file.animation_window_size("user_defined")
  ensight.file.animation_window_xy(1920,1080)
  ensight.file.animation_multiple_images("OFF")
  ensight.file.animation_play_flipbook("OFF")
  ensight.file.animation_play_time("ON")
  ensight.file.animation_reset_flipbook("OFF")
  ensight.file.animation_reset_traces("OFF")
  ensight.file.animation_reset_time("OFF")
  ensight.file.save_animation()
