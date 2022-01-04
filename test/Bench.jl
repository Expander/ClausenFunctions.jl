import ClausenFunctions

n = 10_000_000
x_min = 0.0
x_max = 6.5
data = (x_max - x_min)*rand(Float64, n) + x_min*ones(n)

println("Benchmarking cl2::Float64")

time_cl2(data) = @time map(ClausenFunctions.cl2, data)

map(ClausenFunctions.cl2, data) # trigger compilation
time_cl2(data)                  # trigger compilation
time_cl2(data)
time_cl2(data)

println("Benchmarking cl3::Float64")

time_cl3(data) = @time map(ClausenFunctions.cl3, data)

map(ClausenFunctions.cl3, data) # trigger compilation
time_cl3(data)                  # trigger compilation
time_cl3(data)
time_cl3(data)

println("Benchmarking cl4::Float64")

time_cl4(data) = @time map(ClausenFunctions.cl4, data)

map(ClausenFunctions.cl4, data) # trigger compilation
time_cl4(data)                  # trigger compilation
time_cl4(data)
time_cl4(data)

println("Benchmarking cl5::Float64")

time_cl5(data) = @time map(ClausenFunctions.cl5, data)

map(ClausenFunctions.cl5, data) # trigger compilation
time_cl5(data)                  # trigger compilation
time_cl5(data)
time_cl5(data)

println("Benchmarking cl6::Float64")

time_cl6(data) = @time map(ClausenFunctions.cl6, data)

map(ClausenFunctions.cl6, data) # trigger compilation
time_cl6(data)                  # trigger compilation
time_cl6(data)
time_cl6(data)

println("Benchmarking cl::Float64")

time_cl(k, data) = @time map(x -> ClausenFunctions.cl(k, x), data)

n = 1_000_000
data = (x_max - x_min)*rand(Float64, n) + x_min*ones(n)

for k in [2, 3, 4, 5, 6, 1000, 1001, 1_000_000]
    println("Benchmarking cl($(k),x)::Float64")
    map(x -> ClausenFunctions.cl(k, x), data) # trigger compilation
    time_cl(k, data)            # trigger compilation
    time_cl(k, data)
    time_cl(k, data)
end
