using PyPlot

b_ref = DataFrame(bus = collect(1:15), v_squared = [1.000, 0.942,  0.964,  1.000, 
             0.997,  0.994,  0.992,  1.041, 1.021,  1.023,  
             1.031,  1.034, 0.959,  0.950,  0.944])

function compare_voltages(res_dict1, res_dict2; label1="", label2="" )

u1 = res_dict1[:v_squared]
u2 = res_dict2[:v_squared] 
buses = res_dict1[:bus]

v1 = sqrt.(u1)
v2 = sqrt.(u2)

if length(v1) != length(v1)
    warn("Plotter: length of input vectors not the same")
end

clf()
plot(buses, v1, "ro-", label=label1)
plot(buses, v2, "bo-", label=label2)
legend()

end


function compare_voltages(res_dict1, res_dict2, res_dict3; label1="", label2="", label3="" )

u1 = res_dict1[:v_squared]
u2 = res_dict2[:v_squared]
u3 = res_dict3[:v_squared]
buses = res_dict1[:bus]

v1 = sqrt.(u1)
v2 = sqrt.(u2)
v3 = sqrt.(u3)

if length(v1) != length(v1)
    warn("Plotter: length of input vectors not the same")
end

clf()
plot(buses, v1, "ro-", label=label1)
plot(buses, v2, "bo-", label=label2)
plot(buses, v3, "ko-", label=label3)
legend()

end