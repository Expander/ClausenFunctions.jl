import Clausen

n = 10_000_000
x_min = 0.0
x_max = 6.5
data = (x_max - x_min)*rand(Float64, n) + x_min*ones(n)

println("Benchmarking cl2::Float64")

time_cl2(data) = @time map(Clausen.cl2, data)

map(Clausen.cl2, data)       # trigger compilation
time_cl2(data)               # trigger compilation
time_cl2(data)
time_cl2(data)
