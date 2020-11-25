module KhepriFrame3DD
using KhepriBase

# functions that need specialization
include(khepribase_interface_file())

include("Frame3DD.jl")

function __init__()

#=  set_backend_family(default_wall_family(), autocad, acad_layer_family("Wall"))
  set_backend_family(default_slab_family(), autocad, acad_layer_family("Slab"))
  set_backend_family(default_roof_family(), autocad, acad_layer_family("Roof"))
  set_backend_family(default_beam_family(), autocad, acad_layer_family("Beam"))
  set_backend_family(default_column_family(), autocad, acad_layer_family("Column"))
  set_backend_family(default_door_family(), autocad, acad_layer_family("Door"))
  set_backend_family(default_panel_family(), autocad, acad_layer_family("Panel"))

  set_backend_family(default_table_family(), autocad, acad_layer_family("Table"))
  set_backend_family(default_chair_family(), autocad, acad_layer_family("Chair"))
  set_backend_family(default_table_chair_family(), autocad, acad_layer_family("TableChairs"))

  set_backend_family(default_curtain_wall_family(), autocad, acad_layer_family("CurtainWall"))
  set_backend_family(default_curtain_wall_family().panel, autocad, acad_layer_family("CurtainWall-Panel"))
  set_backend_family(default_curtain_wall_family().boundary_frame, autocad, acad_layer_family("CurtainWall-Boundary"))
  set_backend_family(default_curtain_wall_family().transom_frame, autocad, acad_layer_family("CurtainWall-Transom"))
  set_backend_family(default_curtain_wall_family().mullion_frame, autocad, acad_layer_family("CurtainWall-Mullion"))
  =#
  set_backend_family(default_truss_bar_family(),
     frame3dd,
     frame3dd_truss_bar_family(
       E=210e09,         # E (Young's modulus)
       G=81e09,          # G (Kirchoff's or Shear modulus)
       p=0.0,            # Roll angle
       d=7.701e3))       # RO (Density)

  add_current_backend(frame3dd)
end

end
