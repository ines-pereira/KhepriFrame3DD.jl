export frame3dd,
       frame3DD_circular_tube_truss_bar_family

Base.@kwdef struct Frame3DDBackend{K,T} <: LazyBackend{K,T}
  realized::Parameter{Bool}=Parameter(false)
  truss_nodes::Vector{<:TrussNode}=TrussNode[]
  truss_bars::Vector{<:TrussBar}=TrussBar[]
  shapes::Vector{<:Shape}=Shape[] # This contains all the rest that is not treated yet
  truss_node_data::Vector{TrussNodeData}=TrussNodeData[]
  truss_bar_data::Vector{TrussBarData}=TrussBarData[]
end

abstract type FR3DDKey end
const FR3DDId = Any
const FR3DDIds = Vector{FR3DDId}
const FR3DDRef = GenericRef{FR3DDKey, FR3DDId}
const FR3DDRefs = Vector{FR3DDRef}
const FR3DDEmptyRef = EmptyRef{FR3DDKey, FR3DDId}
const FR3DDUniversalRef = UniversalRef{FR3DDKey, FR3DDId}
const FR3DDNativeRef = NativeRef{FR3DDKey, FR3DDId}
const FR3DDUnionRef = UnionRef{FR3DDKey, FR3DDId}
const FR3DDSubtractionRef = SubtractionRef{FR3DDKey, FR3DDId}
const FR3DD = Frame3DDBackend{FR3DDKey, FR3DDId}

void_ref(b::FR3DD) = FR3DDNativeRef(-1)

const frame3dd = FR3DD()

KhepriBase.backend_name(b::FR3DD) = "Frame3DD"

# Frame3DD needs to merge nodes and bars
save_shape!(b::FR3DD, s::TrussNode) = maybe_merged_node(b, s)
save_shape!(b::FR3DD, s::TrussBar) = maybe_merged_bar(b, s)

# Frame3DD does not need layers
with_family_in_layer(f::Function, backend::FR3DD, family::Family) = f()

# Frame3DD Families
abstract type Frame3DDFamily <: Family end

Base.@kwdef struct Frame3DDTrussBarGeometry
  Ax  # cross section area
  Asy # shear area y-direction
  Asz # shear area z-direction
  Jxx # torsional moment of inertia - x axis
  Iyy # bending moment of inertia - y axis
  Izz # bending moment of inertia - z axis
end

Base.@kwdef struct Frame3DDTrussBarFamily <: Frame3DDFamily
  geometry=Ref{Any}(nothing) # properties that only depend on radius and thickness
  E   # elastic modulus
  G   # shear modulus
  p   # roll angle
  d   # mass density
end

export frame3dd_truss_bar_family
frame3dd_truss_bar_family = Frame3DDTrussBarFamily

truss_bar_family_cross_section_area(f::Frame3DDTrussBarFamily) =
  f.geometry[].Ax

#=
Circular Tube (outer radius= Ro, inner radius = Ri):

Ax = π ( Ro2 - Ri2 )
Asy = Asz = Ax / ( 0.54414 + 2.97294(Ri/Ro) - 1.51899(Ri/Ro)2 ) ± 0.05%
Jxx = (1/2) π ( Ro4 - Ri4 )
Ixx = Iyy = (1/4) π ( Ro4 - Ri4 )
=#

#=
outer radius, thickness;
E   # elastic modulus
G   # shear   modulus
p   # roll angle
d   # mass density
=#
export frame3DD_circular_tube_truss_bar_geometry
frame3dd_circular_tube_truss_bar_geometry(rₒ, e) =
  let rᵢ = rₒ - e,
      Ax = annulus_area(rₒ, rᵢ),
      Asyz = Ax/(0.54414 + 2.97294*(rᵢ/rₒ) - 1.51899*(rᵢ/rₒ)^2),
      Jxx = π/2*(rₒ^4 - rᵢ^4),
      Ixxyy = Jxx/2
    Frame3DDTrussBarGeometry(
      Ax =Ax,
      Asy=Asyz,
      Asz=Asyz,
      Jxx=Jxx,
      Iyy=Ixxyy,
      Izz=Ixxyy)
 end

#=
 @rad1 = 25.0                   # Section diameter (mm)
 @thk1 = 1.6                    # Section thickness (mm)

 @section = "CHS_25x1.6"        # Tube name
 @material = "Steel C450"       # Duragal DualGrade RHS Grade C450Lo (AS 1163)
 @ftype = "Tube"                # Tube type (Tube, Beam, Channel, Angle)
 @emod = "210.0"                # Young's Modulus (GN per square metre) for steel
 @gmod = "79.0"                 # Modulus of Rigidity (GN per square metre) for steel
 @pang = "0.0"                  # section axis rotation angle (deg)
 @ax = "568.3"                  # cross-sectional area (mm^2)
 @asy = "3.2"                   # 2 x thickness (mm) shear area
 @asz = "3.2"                   # 2 x thickness (mm) shear area
 @jxx = "0.0585E+06"            # J (mm^4)
 @iyy = "0.0702E+06"            # Iyy (mm^4)
 @izz = "0.0237E+06"            # Izz (mm^4)
=#

backend_get_family_ref(b::FR3DD, f::TrussBarFamily, tbf::Frame3DDTrussBarFamily) =
  begin
    tbf.geometry[] = frame3dd_circular_tube_truss_bar_geometry(f.radius, f.radius-f.inner_radius)
    tbf
  end

#
KhepriBase.b_delete_all_refs(b::FR3DD) =
  begin
    empty!(b.truss_nodes)
    empty!(b.truss_bars)
    b.realized(false)
    b
  end

realize(b::FR3DD, s::TrussNode) =
  error("BUM")

realize(b::FR3DD, s::TrussBar) =
  error("BUM")

realize(b::FR3DD, s::Panel) =
  error("BUM")


#new_truss_analysis(v=nothing; self_weight=false, backend=frame3dd) =
displacements_from_frame3dd(b::FR3DD, filename, load::Vec, self_weight::Bool, point_loads::Dict) =
    let nodes = process_nodes(b.truss_nodes, load, point_loads),
      bars = process_bars(b.truss_bars, nodes),
      supports = unique(filter(s -> s != false, map(n -> n.family.support, nodes))),
      loaded_nodes = filter(!truss_node_is_supported, nodes),
      supported_nodes = filter(truss_node_is_supported, nodes),
      shear = 0,		     # 1: include shear deformation effects, 0: don't
      geom  = 0,		     # 1: include geometric stiffness effects, 0: don't
      exagg = 1,		     # exaggeration factor for Gnuplot output
      scale = 1.0,		   # zoom scale factor for Gnuplot plotting
      dx = 1.0		       # x-axis increment for internal force calc's
    empty!(b.truss_node_data)
    append!(b.truss_node_data, nodes)
    empty!(b.truss_bar_data)
    append!(b.truss_bar_data, bars)
    open(filename, "w") do io
      println(io, "Frame analysis for Khepri");
      println(io, "%% node data ...");
      println(io, "$(length(nodes))\t\t%% number of nodes");
      println(io, "%% J\t\tX\t\tY\t\tZ\t\tr\t\tnode coordinates");
      for n in nodes
        i = n.id
        p = n.loc
        println(io, "$(i)\t$(p.x)\t$(p.y)\t$(p.z)\t0")
      end
      println(io, "");
      println(io, "%% reaction data ...");
      #nR = count(has_support, nodes);
      println(io, "$(length(supported_nodes))    %% number of nodes with reaction forces");
      println(io, "%% j\tRx\tRy\tRz\tRxx\tRyy\tRzz");
      for n in supported_nodes
        s = n.family.support
        print(io, "$(n.id),\t$(Int(s.ux)),\t$(Int(s.uy)),\t$(Int(s.uz))\t")
        println(io, "$(Int(s.rx)),\t$(Int(s.ry)),\t$(Int(s.rz))")
      end
      println(io, "")
      println(io, "%% frame element section property data ...")
      println(io, "$(length(bars))\t\t%% number of frame elements")
      println(io, "%% m\tn1\tn2\t\tAx\t\tAsy\t\tAsz\t\tJxx\t\tIyy\t\tIzz\t\tE\t\tG\t\tp\tdensity");
      for bar in bars
        f = family_ref(b, bar.family)
        g = f.geometry[]
        print(io, "$(bar.id),\t$(bar.node1.id),\t$(bar.node2.id),")
        print(io, "\t$(g.Ax),\t$(g.Asy),\t$(g.Asz),")
        print(io, "\t$(g.Jxx),\t$(g.Iyy),\t$(g.Izz),")
        println(io, "\t$(f.E),\t$(f.G),\t$(f.p),\t$(f.d)")
      end
      println(io, "\n");
      println(io, "$(shear)\t\t%% 1: include shear deformation, 0: do not");
      println(io, "$(geom)\t\t%% 1: include geometric stiffness, 0: do not");
      println(io, "$(exagg)\t%% exagerate deformations in plotting");
      println(io, "$(scale)\t%% zoom scale factor for 3D plotting");
      println(io, "$(dx)\t%% x-axis increment for internal forces calc");
      println(io, "\n");
      println(io, "%% static load data ...");
      println(io, "1\t\t%% number of static load cases");
      println(io, "\t\t%% begin static load case 1 of 1");
      println(io, "%% gravitational acceleration for self-weight loading");
      println(io, "%% gX         gY         gZ");
      println(io, "  0.0        0.0        $(self_weight ? -9.81 : 0.0)\n");
      println(io, "$(length(loaded_nodes))\t\t%% number of loaded nodes");
      println(io, "%% j\t\tFx\t\tFy\t\tFz\t\tMxx\t\tMyy\t\tMzz");
      for n in loaded_nodes
        println(io, "$(n.id),\t$(n.load.x),\t$(n.load.y),\t$(n.load.z),\t0,\t0,\t0")
      end
      println(io, "\n");
      println(io, "0\t\t%% number of members with uniform distributed loads")
      println(io, "0\t\t%% number of members with trapezoidal loads");
      println(io, "0\t\t%% number of members with internal point loads");
      println(io, "0\t\t%% number of members with temperature loads");
      println(io, "0\t\t%% number of nodes with prescribed displacements");
      println(io, "\t\t%% end   static load case 1 of 1\n");
      println(io, "%% inertial load data ...");
      println(io, "0\t\t%% number of dynamic modes to analyze");
    end
  end

frame3dd_plugin = joinpath(dirname(abspath(@__DIR__)), "bin", "frame3dd.exe")

frame3dd_simulation_path() =
  mktempdir(tempdir(), prefix="Frame3DD_")


KhepriBase.b_truss_analysis(b::FR3DD, load::Vec, self_weight::Bool, point_loads::Dict) =
  # Ensure extension is FMM to force Matlab mode
  let simulation_folder = frame3dd_simulation_path(),
       input_path = joinpath(simulation_folder, "IOdata.IN"),
       output_path = joinpath(simulation_folder, "IOdata.OUT")
    #@info input_path
    displacements_from_frame3dd(b, input_path, load, self_weight, point_loads)
    #try
      withenv("FRAME3DD_OUTDIR"=>simulation_folder) do
        run(`$frame3dd_plugin -q -i $input_path -o $output_path`)
      end
    #catch e
    #end
    lines = readlines(output_path)
    idx1 = findfirst(l->startswith(l, "N O D E   D I S P L A C E M E N T S"), lines) + 2
    idx2 = findnext(l->startswith(l, "F R A M E   E L E M E N T   E N D   F O R C E S"), lines, idx1)
    Dict([let strs = split(line)
            parse(Int, strs[1])=>parse.(Float64, strs[2:7])
          end
          for line in lines[idx1:idx2-1]])
  end

KhepriBase.b_node_displacement_function(b::FR3DD, results) =
  n -> vxyz(get(results, n.id, (0,0,0,0,0,0))[1:3]..., world_cs)

#=
using Libdl
f3ddlib = dlopen(joinpath(dirname(abspath(@__DIR__)), "bin", "Frame3DDLib"))
f3ddmain = dlsym(f3ddlib, :simulate)
=#

#dlclose(f3ddlib)

const f3ddlibpath = joinpath(dirname(abspath(@__DIR__)), "bin", "Frame3DDLib")
##########################################
KhepriBase.b_truss_analysis(be::FR3DD, load::Vec, self_weight::Bool, point_loads::Dict) =
    let nodes = process_nodes(be.truss_nodes, load, point_loads),
        bars = process_bars(be.truss_bars, nodes),
        supports = unique(filter(s -> s != false, map(n -> n.family.support, nodes))),
        loaded_nodes = filter(!truss_node_is_supported, nodes),
        supported_nodes = filter(truss_node_is_supported, nodes),
        shear = 0,		     # 1: include shear deformation effects, 0: don't
        geom  = 0,		     # 1: include geometric stiffness effects, 0: don't
        exagg = 1,		     # exaggeration factor for Gnuplot output
        scale = 1.0,		   # zoom scale factor for Gnuplot plotting
        dx = 1.0,		       # x-axis increment for internal force calc's
        nN = length(nodes), # number of nodes
        nE = length(bars),  # number of bars
        nodes_coords = zeros(Float64, 3*(nN+1)), # X,Y,Z node coordinates (global)
        nodes_radius = zeros(Float32, nN + 1),   # node size radius, for finite sizes
        DoF = 6*nN,
        nR = length(supported_nodes),            # number of restrained nodes
        q = zeros(Int32, DoF+1),
        r = zeros(Int32, DoF+1),
        sumR = 0
        empty!(be.truss_node_data)
        append!(be.truss_node_data, nodes)
        empty!(be.truss_bar_data)
        append!(be.truss_bar_data, bars)
      for n in nodes
          i = n.id
          p = n.loc
        for j in 1:3
          nodes_coords[3*i+j]=p.raw[j]
        end
      end
      for n in supported_nodes
        s = n.family.support
        j = n.id
        sj = [Int(s.ux), Int(s.uy), Int(s.uz), Int(s.rx), Int(s.ry), Int(s.rz)]
        for l in 5:-1:0
          r[6*j-l+1] = sj[6-l]
        end
      end
      sumR = sum(r)
      sumR > 4 || error("Un-restrained structure: only $(sumR) imposed reactions. At least 4 reactions are required to support static loads.")
      sumR < DoF || error("Fully restrained structure: $(sumR) imposed reactions >= $(DoF) degrees of freedom")
      for i in 1:DoF
        q[i+1] = r[i+1] > 0 ? 0 : 1
      end
      L = zeros(Float64, nE+1)   # length of each element
      Le = zeros(Float64, nE+1)  # effective length of each element
      N1  = zeros(Int32, nE+1)   # node #1 of each element
      N2  = zeros(Int32, nE+1)   # node #2 of each element
    	Ax  = zeros(Float32, nE+1) # cross section area of each element	*/
    	Asy = zeros(Float32, nE+1) # shear area in local y direction 	*/
    	Asz = zeros(Float32, nE+1) # shear area in local z direction	*/
    	Jx  = zeros(Float32, nE+1) # torsional moment of inertia 		*/
    	Iy  = zeros(Float32, nE+1) # bending moment of inertia about y-axis */
    	Iz  = zeros(Float32, nE+1) # bending moment of inertia about z-axis */
    	E   = zeros(Float32, nE+1) # frame element Young's modulus	*/
    	G   = zeros(Float32, nE+1) # frame element shear modulus		*/
    	p   = zeros(Float32, nE+1) # element rotation angle about local x axis */
    	d   = zeros(Float32, nE+1) # element mass density			*/
      for bar in bars
        f = family_ref(be, bar.family)
        g = f.geometry[]
        b = bar.id
        N1[b+1] = bar.node1.id
        N2[b+1] = bar.node2.id
        Ax[b+1] = g.Ax
        Asy[b+1] = g.Asy
        Asz[b+1] = g.Asz
        Jx[b+1] = g.Jxx
        Iy[b+1] = g.Iyy
        Iz[b+1] = g.Izz
        E[b+1] = f.E
        G[b+1] = f.G
        p[b+1] = f.p
        d[b+1] = f.d
        n1 = N1[b+1]
        n2 = N2[b+1]
        L[b+1] = distance(bar.node2.loc, bar.node1.loc)
        Le[b+1] = L[b+1] - nodes_radius[n2+1] - nodes_radius[n1+1]
      end
      nF = length(loaded_nodes)
      F_mech = zeros(Float64, DoF+1)
      for n in loaded_nodes
        F_mech[6*n.id - 5 + 1] = n.load.x;
        F_mech[6*n.id - 4 + 1] = n.load.y;
        F_mech[6*n.id - 3 + 1] = n.load.z;
      end
      gX = [0f0, 0f0]
      gY = [0f0, 0f0]
      gZ = [0f0, self_weight ? -9.81f0 : 0f0]
      nFs = Int32[0, nF]
      F_mechs = [Ref(F_mech, 1), Ref(F_mech, 1)]
      D = zeros(DoF+1)	  # displacements of each node
      #ccall(f3ddmain, Int32, (
      #      Int32,        # number of Nodes
      #      Int32,        # number of frame Elements
      #      Ref{Float64}, # X,Y,Z node coordinates (global)
      #      Ref{Float32}, # node size radius, for finite sizes
      #      Int32,        # number of restrained nodes
      #      Ref{Int32},   # reaction data
      #      Ref{Int32},   # reaction data
      #      Int32,        # total no. of reactions
      #      Ref{Float64},	# length of each element
      #      Ref{Float64},	# effective length of each element
      #      Ref{Int32},   # node #1 of each element
      #      Ref{Int32},   # N2,    // node #2 of each element
      #      Ref{Float32},	# cross section area of each element
      #      Ref{Float32},	# shear area in local y direction 	*/
      #      Ref{Float32},	# shear area in local z direction	*/
      #      Ref{Float32},	# torsional moment of inertia 		*/
      #      Ref{Float32},	# bending moment of inertia about y-axis */
      #      Ref{Float32},	# bending moment of inertia about z-axis */
      #      Ref{Float32},	# frame element Young's modulus	*/
      #      Ref{Float32},	# frame element shear modulus		*/
      #      Ref{Float32},	# element rotation angle about local x axis */
      #      Ref{Float32},	# element mass density			*/
      #      Int32,        # number of Load cases
      #      Int32,        # number of Load cases
      #      Ref{Float32}, # gravitational acceleration in global X
      #      Ref{Float32}, # gravitational acceleration in global Y
      #      Ref{Float32}, # gravitational acceleration in global Z
      #      Ref{Int32},   # number of loaded nodes
      #      Ref{Ptr{Float64}},  # mechanical load vectors, all load cases
      #      # RESULTS
      #      Ref{Float64}  # displacements of each node
      #      ),
      #  nN, nE,
      #  nodes_coords, nodes_radius,
      #  nR, q, r, sumR,
      #  L, Le, N1, N2, Ax, Asy, Asz, Jx, Iy, Iz, E, G, p, d,
      #  1, 0,
      #  gX, gY, gZ,
      #  nFs,
      #  F_mechs,
      #  # RESULTS
      #  D)
      @ccall f3ddlibpath.simulate(
         nN::Int32,        # number of Nodes
         nE::Int32,        # number of frame Elements
         nodes_coords::Ref{Float64}, # X,Y,Z node coordinates (global)
         nodes_radius::Ref{Float32}, # node size radius, for finite sizes
         nR::Int32,        # number of restrained nodes
         q::Ref{Int32},   # reaction data
         r::Ref{Int32},   # reaction data
         sumR::Int32,        # total no. of reactions
         L::Ref{Float64},	# length of each element
         Le::Ref{Float64},	# effective length of each element
         N1::Ref{Int32},   # node #1 of each element
         N2::Ref{Int32},   # node #2 of each element
         Ax::Ref{Float32},	# cross section area of each element
         Asy::Ref{Float32},	# shear area in local y direction 	*/
         Asz::Ref{Float32},	# shear area in local z direction	*/
         Jx::Ref{Float32},	# torsional moment of inertia 		*/
         Iy::Ref{Float32},	# bending moment of inertia about y-axis */
         Iz::Ref{Float32},	# bending moment of inertia about z-axis */
         E::Ref{Float32},	# frame element Young's modulus	*/
         G::Ref{Float32},	# frame element shear modulus		*/
         p::Ref{Float32},	# element rotation angle about local x axis */
         d::Ref{Float32},	# element mass density			*/
         1::Int32,        # number of Load cases
         0::Int32,        # number of Load cases
         gX::Ref{Float32}, # gravitational acceleration in global X
         gY::Ref{Float32}, # gravitational acceleration in global Y
         gZ::Ref{Float32}, # gravitational acceleration in global Z
         nFs::Ref{Int32},   # number of loaded nodes
         F_mechs::Ref{Ptr{Float64}},  # mechanical load vectors, all load cases
         # RESULTS
         D::Ref{Float64}  # displacements of each node
      )::Int32
      D
    end
#  end

KhepriBase.b_node_displacement_function(b::FR3DD, results) =
  n -> vxyz(results[6*n.id-4:6*n.id-2]..., world_cs)
