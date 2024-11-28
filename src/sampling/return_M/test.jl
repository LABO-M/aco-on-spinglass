using ProgressMeter

n = 10
progress = Progress(n)

for i in 1:n
    next!(progress)
    sleep(0.1)
    next!(progress)
end