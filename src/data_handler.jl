using DataFrames

# TYPE DEFINITIONS
type Generator
   index::Any
   bus_idx::Int
   g_P_max::Float64
   g_Q_max::Float64
   cost::Float64
   function Generator(index, bus_idx, g_P_max, g_Q_max)
      g = new()
      g.index  = index
      g.bus_idx = bus_idx
      g.cost = 1.00
      g.g_P_max = g_P_max
      g.g_Q_max = g_Q_max
      return g
   end
end

type Bus
   index::Any
   is_root::Bool
   d_P::Float64
   d_Q::Float64
   cosphi::Float64
   tanphi::Float64
   v_max::Float64
   v_min::Float64
   children::Vector{Int}
   ancestor::Vector{Int}
   generator::Generator
   function Bus(index, d_P, d_Q, v_max, v_min)
      b = new()
      b.index = index
      b.is_root = false
      b.d_P = d_P
      b.d_Q = d_Q
      b.v_max = v_max
      b.v_min = v_min
      b.children = Int[]
      b.ancestor = Int[]
      cosphi = d_P/(sqrt(d_P^2 + d_Q^2))
      tanphi = tan(acos(b.cosphi)) 
      if isnan(cosphi)
        b.cosphi = 0
        b.tanphi = 0
      else
        b.cosphi = cosphi
        b.tanphi = tan(acos(cosphi)) 
      end
      return b
   end
end

type Line
   index::Any
   to_node::Int # the "to" node
   from_node::Int # the "from" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   b::Float64 # the susceptance value
   s_max::Float64 # the capacity of the line
   function Line(index, to_node, from_node, r, x, s_max)
      l = new()
      l.index = index
      l.to_node = to_node
      l.from_node = from_node
      l.r = r
      l.x = x
      l.b = (x/(r^2 + x^2))
      l.s_max = s_max
      return l
   end
end


function load_case_data(;datafile="")
# READ RAW DATA
info("Reading Data")

if datafile == ""
  data_dir = DATA_DIR
else
  data_dir = datafile
end

nodes_raw = readtable("data/$data_dir/nodes.csv")
sum(nonunique(nodes_raw, :index)) != 0 ? warn("Ambiguous Node Indices") : nothing

lines_raw = readtable("data/$data_dir/lines.csv")
sum(nonunique(lines_raw, :index)) != 0  ? warn("Ambiguous Line Indices") : nothing

generators_raw = readtable("data/$data_dir/generators.csv")
sum(nonunique(generators_raw, :index)) != 0 ? warn("Ambiguous Generator Indices") : nothing


# PREPARING MODEL DATA

buses = Dict()
for n in 1:nrow(nodes_raw)
    index = nodes_raw[n, :index]
    d_P = nodes_raw[n, :d_P]
    d_Q = nodes_raw[n, :d_Q]
    v_max = nodes_raw[n, :v_max]
    v_min = nodes_raw[n, :v_min]
    newb = Bus(index, d_P, d_Q, v_max, v_min)
    buses[newb.index] = newb
end

lines = Dict()
for l in 1:nrow(lines_raw)
    index = lines_raw[l, :index]
    from_node = lines_raw[l, :from_node]
    to_node = lines_raw[l, :to_node]
    r = lines_raw[l, :r]
    x = lines_raw[l, :x]
    b = lines_raw[l, :b]
    s_max = lines_raw[l, :s_max]
    newl = Line(index, to_node, from_node, r, x, s_max)

    push!(buses[newl.from_node].children, newl.to_node)
    push!(buses[newl.to_node].ancestor, newl.from_node)

    lines[newl.index] = newl
end

# Check topology
r = 0
root_bus = 0
for b in keys(buses)
    l = length(buses[b].ancestor)
    if l > 1
        warn("Network not Radial (Bus $(buses[b].index))")
    elseif l == 0
        buses[b].is_root = true
        root_bus = b
        r += 1
    end
end
if r == 0
    warn("No root detected")
    root_bus = 0
elseif r > 1
    warn("More than one root detected")
end

generators = Dict()
for g in 1:nrow(generators_raw)
    index = generators_raw[g, :index]
    bus_idx = generators_raw[g, :node]
    g_P_max = generators_raw[g, :p_max]
    g_Q_max = generators_raw[g, :q_max]
    cost = generators_raw[g, :cost]
    newg = Generator(index, bus_idx, g_P_max, g_Q_max)
    newg.cost = cost

    buses[newg.bus_idx].generator = newg

    generators[newg.index] = newg
end

# Radial PTDF
A = zeros(length(lines), length(buses))
for b in collect(keys(buses))
    a = b
    while a != root_bus
        A[a-1, b] = 1
        a = buses[a].ancestor[1]
    end
end

info("Done preparing Data")
return buses, lines, generators, A
end


function change_load_same_pf(α)
# always uses the GLOBAL BUSES as reference

  buses = deepcopy(BUSES)
  for i in keys(buses)
      P = buses[i].d_P
      Q = buses[i].d_Q
      PF = P/sqrt(P^2 + Q^2)
      if !(isnan(PF))
         buses[i].d_P = P*α
         buses[i].d_Q = P*α*sqrt((1-PF^2)/PF^2)
      end
   end

   return buses
end


function add_generator(index, bus, g_P_max, g_Q_max, cost)
    newg = Generator(index, bus, g_P_max, g_Q_max)
    newg.cost = cost
    BUSES[newg.bus_idx].generator = newg
    GENERATORS[newg.index] = newg
end
