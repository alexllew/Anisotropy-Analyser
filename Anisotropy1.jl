
using ExcelReaders, DelimitedFiles, DSP, Plots, Statistics, LinearAlgebra,  Spectra, Dierckx, Glob

fileNames = glob("*.xlsx", "Data")

const min = 0.1
const max = 100
const delta = 0.0195
const M = 1000
const N = 3000
const eta = 2
const sigma = 0.1
const L = 10
const alpha_max = 10^30
const alpha_min = 10^-30
const alpha_0 = 1
const start_frac = 0.25
const pen = 1
const tolP = 10^-8
const max_iterations = 10000
const report = 1000
const lag = 5
const threshold = 5
const influence = 0
const backrange = 200
const minback = 10^-10

function optimise()
    for a in 1:length(fileNames)
        println("Initialising...")
        fileName="$(fileNames[a])"
		if isdir(splitext(fileName)[1]) != true
			mkdir(splitext(fileName)[1])
		end
		A, A_r, Areas, Areas_r, Sm, Dm, bS, bD, f_prev, Taus, Time, logyfacS, logyfacD = initialise(a, fileName)
		solution_history = []
		Params = A, Sm, bS
		solution, model = Poisson(f_prev, Params)
		push!(solution_history, solution)
		f_now = next_iter(f_prev, alpha_0, model, Params)
		solution, model = Poisson(f_now, Params)
		push!(solution_history, solution)
		escape = 0
		criterion = 0
		i = 0
		for i in 3:max_iterations
			d_k = f_now - f_prev
			alpha_k = alpha_init(d_k, model, Params)
			f_next = next_iter(f_now, alpha_k, model, Params)
			solution, poss_model = Poisson(f_next, Params)
			alpha_count = 0
			while acceptance(solution, solution_history, alpha_k, f_next, f_now) == 0
				if alpha_count < 100
					alpha_k = incremalpha(alpha_k)
					f_next = next_iter(f_now, alpha_k, model, Params)
					solution, poss_model = Poisson(f_next, Params)
					alpha_count+=1
				else
					escape = 1
					break
				end
			end
			if escape == 1
				break
			end
			model = poss_model
			if length(solution_history)>L
				popfirst!(solution_history)
			end
			push!(solution_history, solution)
			if mod(i, report) == 0
				criterion = check_termination(f_next, f_now)
				if criterion < tolP
					f_now = f_next
					break
				end
				tracedata(fileName ,solution, criterion, logyfacS, i, Time, model, Taus, f_now, Sm, bS, Areas, 0)
			end
			f_prev = f_now
			f_now = f_next
		end
		tracedata(fileName ,solution, criterion, logyfacS, i, Time, model, Taus, f_now, Sm, bS, Areas, 0, 0, 1)
		f_lt = f_now
		lifetimeweight = f_lt.*Areas
		writedlm("$(splitext(fileName)[1])/Lifetime.csv", [f_lt, lifetimeweight], ",")
		println("Done")
		optimiseDm(f_lt, Dm, Taus, bD, A, A_r, Areas, Areas_r, logyfacD, Time, fileName)
    end
end

function optimiseDm(f_lt, Dm, Taus, bD, A_old, A, Areas, Areas_r, logyfac, Time, fileName)
    f_prev, oldmodel = f_init(A_old, A, f_lt, Dm)
    solution_history = []
    Params = A, Dm, bD
    solution, model = Poisson(f_prev, oldmodel, Params)
    push!(solution_history, solution)
    f_now = next_iter(f_prev, alpha_0, model, oldmodel, Params)
    solution, model = Poisson(f_now, oldmodel, Params)
    push!(solution_history, solution)
    criterion = 1
    escape = 0
    i=0
    for i in 3:max_iterations*100
        d_k = f_now - f_prev
        alpha_k = alpha_init(d_k, oldmodel, model, Params)
        f_next = next_iter(f_now, alpha_k, model, oldmodel, Params)
        solution, poss_model = Poisson(f_next, oldmodel, Params)
        alpha_count = 0
        while acceptance(solution, solution_history, alpha_k, f_next, f_now) == 0
            if alpha_count<100
                alpha_k = incremalpha(alpha_k)
                f_next = next_iter(f_now, alpha_k, model, oldmodel, Params)
                solution, poss_model = Poisson(f_next, oldmodel, Params)
                alpha_count+=1
            else
                escape = 1
                break
            end
        end
        if escape == 1
            break
        end
        model = poss_model
        if length(solution_history)>L
            popfirst!(solution_history)
        end
        push!(solution_history, solution)
        if mod(i, report) == 0
            criterion = check_termination(f_next, f_now)
            if criterion < tolP
                f_now = f_next
                tracedata(fileName, solution, criterion, logyfac, i, Time, model, Taus, f_now, Dm, bD, Areas_r, 1, A)
                break
            end
            tracedata(fileName, solution, criterion, logyfac, i, Time, model, Taus, f_next, Dm, bD, Areas_r, 1, A)
        end
        f_prev = f_now
        f_now = f_next
    end
    f_r_final = f_now
    tracedata(fileName, solution, criterion, logyfac, i, Time, model, Taus, f_r_final, Dm, bD, Areas_r, 1, A, 1)
    lifetimedist = f_lt
    lifetimeweight = f_lt.*Areas
    anisodist = f_r_final
    anisoweight = f_r_final.*Areas_r
    sumdecay = A_old*f_lt
    anisodecay = A*f_r_final
    diffdecay = sumdecay.*anisodecay
    writedlm("$(splitext(fileName)[1])/Lifetime.csv", [lifetimedist, lifetimeweight], ",")
    writedlm("$(splitext(fileName)[1])/Anisotropy.csv", [anisodist, anisoweight], ",")
    writedlm("$(splitext(fileName)[1])/Anisotropy Decay.csv", anisodecay, ",")
    writedlm("$(splitext(fileName)[1])/Decays.csv", [sumdecay, diffdecay], ",")
end

function Poisson(f, Params)
    A, Sm, bS = Params
    bS = bS[1]
    total = 0
    model = A*f
    for i in 1:length(model)
        total += (Sm[i]*log(model[i].+bS))
    end
    solution = sum(model.+bS) - total + pen*sum(f)
    solution, model
end

function Poisson(f_r, oldmodel, Params)
    A, Dm, bD = Params
    bD = bD[1]
    total = 0
    model = oldmodel.*(A*f_r)
    for i in 1:length(model)
        total += (Dm[i]*log(model[i].+bD))
    end
    solution = sum(model.+bD) - total + pen*sum(f_r)
    solution, model
end

function next_iter(f_now, alpha_k, model, Params)
    A, y, b = Params
    grad = Gradient(model, Params)
    f_next = f_now - (1.0/alpha_k) * grad .- pen/alpha_k
    f_next = [if i < 0 0.0 else i end for i in f_next]
    f_next
end

function next_iter(f_now, alpha_k, model, oldmodel, Params)
    A, y, b = Params
    grad = Gradient(model, oldmodel, Params)
    f_next = f_now - (1.0/alpha_k) * grad .- pen/alpha_k
    last = f_next[end]
    f_next = [if i < 0 0.0 else i end for i in f_next]
    f_next
end

function Gradient(model,Params)
    A, y, b = Params
    total = zeros(M)
    for i in 1:length(model)
        total += (A[i,:]).*(ones(M)-((y[i]/(model[i].+b))*ones(M)))
    end
    total
end

function Gradient(model, oldmodel, Params)
    A, y, b = Params
    total = zeros(M)
    for i in 1:length(model)
        total += (oldmodel[i].*A[i,:]).*(ones(M)-((y[i]/(model[i].+b))*ones(M)))
    end
    total
end

function alpha_init(d_k, model, Params)
    A, y, b = Params
    alpha = (norm(((sqrt.(y).*(A*d_k))./(model.+b)), 2)^2)/
        (norm(d_k, 2)^2)
    if alpha < alpha_min
        alpha = alpha_min
    elseif alpha > alpha_max
        alpha = alpha_max
    end
    alpha
end

function alpha_init(d_k, oldmodel, model, Params)
	A, y, b = Params
    alpha = (norm(((sqrt.(y).*oldmodel.*(A*d_k))./(model.+b)), 2)^2)/
        (norm(d_k, 2)^2)
    if alpha < alpha_min
        alpha = alpha_min
    elseif alpha > alpha_max
        alpha = alpha_max
    end
    alpha
end

function acceptance(solution, solution_history, alpha_k, f_next, f_now)
    if solution <= maximum(solution_history)-(((sigma*alpha_k)/2)*
        (norm(f_next - f_now)^2))
        1
    else
        0
    end
end

function incremalpha(alpha_k)
    alpha_k = eta * alpha_k
    if alpha_k < alpha_min
        alpha_k = alpha_min
    elseif alpha_k > alpha_max
        alpha_k = alpha_max
    end
    alpha_k
end

function check_termination(f_next, f_now)
    criterion = norm((f_next.-f_now),2)/norm(f_now, 2)
end

#Data Initialisation

function initialise(a, fileName)
    Taus = Taus_init()
    A, A_r = A_init(fileName, Taus)
    Sm, Dm = y_init(a, fileName)
    f = f_init(A, Sm)
    bS, bD = b_init(Sm, Dm)
    Time = collect(0:delta:delta*(N-1))
    A, A_r, Sm, Dm, Time = reducerange(A, A_r, Sm, Dm, Time)
    Areas, Areas_r = Areas_init(A, A_r, Time)
    logyfacS = Logyfac(Sm)
    logyfacD = Logyfac(Dm)
    A, A_r, Areas, Areas_r, Sm, Dm, bS, bD, f, Taus, Time, logyfacS, logyfacD
end

#Background of y
function b_init(Sm, Dm)
    bS = mean(Sm[1:backrange])
    #bD = mean(Dm[1:backrange])
    bD=1
    bS, bD
end

#Generates logarithmic array of lifetime components from min to max
function Taus_init()
    ratio = (max/min)^(1/(M-1))
    taus = []
    for i in 0:(M-1)
        push!(taus, min*ratio^i)
    end
    taus
end

function y_init(a, fileName)
    ypar = readxl(fileName, "Sheet1!B1:B$N")
    yperp = readxl(fileName, "Sheet1!C1:C$N")
    Sm = ypar .+ 2*yperp
    Dm = ypar .- yperp
    lowest = minimum(Dm)
    if lowest < 0
        Dm = Dm .- lowest .+ 1
    end
    Dm = [if i < 1 1 else i end for i in Dm]
    Sm, Dm
end

#Opens the IRF and converts to a Float64 n-vector
function H(fileName)
	print(fileName)
    h = readxl(fileName, "Sheet1!A1:A$N")
    h = vec(h)
    h = convert(Array{Float64,1}, h)
end

#Generates matrix A, where each column is the convolution of the IRF with
#the exponential decay function to each lifetime component of taus
function A_init(fileName, Taus)
    h = H(fileName)
    decay = []
    Conv = []
    for i in 1:N
        push!(decay, exp((-i*delta)/min))
        decay = convert(Array{Float64,1}, decay)
    end
    A = conv(h, decay)
    A_r = decay
    A = A[1:N,:]
    A_r = A_r[1:N,:]
    for j in 2:(M-1)
        decay = []
        Conv = []
        for i in 1:N
            push!(decay, exp((-i*delta)/Taus[j]))
            decay = convert(Array{Float64,1}, decay)
        end
        Conv = conv(h, decay)[1:N]
        A = hcat(A, Conv)
        A_r = hcat(A_r, decay)

    end
    A = hcat(A, ones(N))
    A_r = hcat(A_r, ones(N))
    A, A_r
end

function Logyfac(y)
    total = 0
    for i in 1:length(y)
        logyfac = (y[i] + 1/2)*log(abs(y[i])) - y[i]  + log(2*π)/2 + 1/(12*y[i])
            total += logyfac
    end
    total
end

function f_init(A, y)
    f = ones(M)
    y_guess = A*f
    f = f*(maximum(y)/maximum(y_guess))
end

function f_init(A_old, A, f_lt, Dm)
    f_r = ones(M)
    oldmodel = A_old*f_lt
    Dm_guess = oldmodel.*(A*f_r)
    f_r = f_r*(maximum(Dm)/maximum(Dm_guess))
    f_r, oldmodel
end

function reducerange(A, A_r,Sm, Dm, Time)
    startvalue = start_frac*maximum(Sm)
    if startvalue == 0
        first = default_first
    else
        first = findfirst(x->x>startvalue, Sm)[1]
    end
    A = A[first:N, :]
    A_r = A_r[first:N, :]
    Sm = Sm[first:N]
    Dm = Dm[first:N]
    Time = Time[first:N]
    A, A_r, Sm, Dm, Time
end

function Areas_init(A, A_r, Time)
    Areas = []
    Areas_r = []
    for i in 1:M
        Area = trapz(Time, A[:,i])
        Area_r = trapz(Time, A_r[:,i])
        push!(Areas, Area)
        push!(Areas_r, Area_r)
    end
    Areas, Areas_r
end

#Print Results

function tracedata(fileName, solution, criterion, logyfac, i, Time, model, Taus, f_now, y, b, Areas, Type, A=0, write=0)
    loglike = solution + logyfac - pen*sum(f_now)
    println("Iteration no: $i, Termination Criterion = $(round(criterion, sigdigits=4))")
    println("Poisson Log-Likelihood = $(round(loglike, sigdigits = 5))")
    PlotData(Type, Time, model, b, Taus, f_now, y, A)
    PrintPeaks(Taus, f_now, Areas, write, A, fileName)
end

function PlotData(Type, Time, model, b, Taus, f_now, y, A)
    if Type == 0
        Title = "Lifetime Distribution"
        Xlabel = "Lifetime"
    else
        Title = "Anisotropy Distribution"
        Xlabel = "Anistropy"
    end
    p1 = plot(Time, model.+b, xlabel = "Time/ns",
        ylabel = "Counts", title = "Lifetime Decay")
    plot!(p1, Time, y)
    if Type == 1
        p2 = plot(Time, A*f_now .- f_now[end],
            xlabel = "Time/ns", ylabel = "Anisotropy",
            title = "Anisotropy Decay")
        plt = plot(p1,p2, layout=(2,1),legend=false)
    else
        plt = p1
    end
    display(plt)
end

#Peak Analysis

function PrintPeaks(Taus, f_now, Areas, write, A, fileName)
    signals = FindPeaks(Taus, f_now)
    numpeaks, peaks, indices = QuantifyPeaks(signals, f_now)
    if length(indices)>0
        if indices[1]==0
            indices[1]=1
        end
    end
    println("Number of peaks: $numpeaks")
    if numpeaks == 0
        return 1
    end
    peakdata=zeros(3, numpeaks)
    yield = zeros(numpeaks)
    for i in 1:numpeaks
        println("Peak $i")
        peak = peaks[i]
        x = convert(Array{Float64,1}, Taus[indices[i]:indices[i]+length(peak)-1])
        peak = convert(Array{Float64,1}, peak)
        xval = 0
        if length(peak) > 3
            spl = Spline1D(x, peak)
            xval = round((x[1]+x[end])/2, digits = 3)
            yval = 0
            found = 0
            while found == 0
                yval =  spl(xval)
                if spl(xval+0.001)>yval
                    xval += 0.001
                elseif spl(xval-0.001)>yval
                    xval -= 0.001
                else
                    found = 1
                end
            end
            highest = xval
            xval = round(x[1], digits = 3)
            while found == 1
                if spl(xval+0.001)<yval/2
                    xval += 0.001
                else
                    found = 2
                end
            end
            startvalue = xval
            xval = round(x[end], digits = 3)
            while found == 2
                if spl(xval-0.001)<yval/2
                    xval -= 0.001
                else
                    found = 3
                end
            end
            endvalue = xval
            fwhm = endvalue-startvalue
        else
            highest = x[2]
            if length(x)>2
                fwhm = (x[3]-x[1])/2
            else
                fwhm = 0
            end
        end
        println("Maximum = $(round(highest, digits = 3)) ns")
        println("FWHM = $(round(fwhm, digits = 3)) ns")
        peakdata[1, i]=highest
        peakdata[2, i]=fwhm
        for j in 2:length(peak)-1
            yield[i] += peak[j]*Areas[indices[i]+j-1]
        end
    end
    totalyield = sum(yield)
    for i in 1:numpeaks
        percentage = (yield[i]/totalyield)*100
        println("Peak $i Contribution: $(round(percentage, sigdigits = 4))%")
        peakdata[3, i]=percentage
    end
    if write == 1
        if A == 0
			println("$fileName")
            writedlm("$(splitext(fileName)[1])/peaks_lifetime.csv", peakdata, ",")
        else
            writedlm("$(splitext(fileName)[1])/peaks_anisotropy.csv", peakdata, ",")
        end
    end
end

function FindPeaks(x,y)
    results = SmoothedZscoreAlgo(y, lag, threshold, influence)
    upper_bound = results[:avgFilter] + threshold * results[:stdFilter]
    lower_bound = results[:avgFilter] - threshold * results[:stdFilter]
    yplot = plot(x,y,color="blue", label="Y",legend=:topleft)
    yplot = plot!(x,upper_bound, color="green", label="Upper Bound",legend=:topleft)
    yplot = plot!(x,results[:avgFilter], color="cyan", label="Average Filter",legend=:topleft)
    yplot = plot!(x,lower_bound, color="green", label="Lower Bound",legend=:topleft)
    signalplot = plot(x,results[:signals],color="red",label="Signals",legend=:topleft)
    plt = plot(yplot,signalplot,layout=(2,1),legend=:topleft)
    display(plt)
    results[:signals]
end

function SmoothedZscoreAlgo(y, lag, threshold, influence)
    padding = zeros(100)
    y=vcat(padding, y)
    n = length(y)
    signals = zeros(n) # init signal results
    filteredY = copy(y) # init filtered series
    avgFilter = zeros(n) # init average filter
    stdFilter = zeros(n) # init std filter
    avgFilter[lag - 1] = mean(y[1:lag]) # init first value
    stdFilter[lag - 1] = std(y[1:lag]) # init first value

    for i in range(lag, stop=n-1)
        if abs(y[i] - avgFilter[i-1]) > threshold*stdFilter[i-1]
            if y[i] > avgFilter[i-1]
                signals[i] += 1 # postive signal
            else
                signals[i] += -1 # negative signal
            end
            # Make influence lower
            filteredY[i] = influence*y[i] + (1-influence)*filteredY[i-1]
        else
            signals[i] = 0
            filteredY[i] = y[i]
        end
        avgFilter[i] = mean(filteredY[i-lag+1:i])
        stdFilter[i] = std(filteredY[i-lag+1:i])
    end
    signals=signals[101:end]
    avgFilter=avgFilter[101:end]
    stdFilter=stdFilter[101:end]
    return (signals = signals, avgFilter = avgFilter, stdFilter = stdFilter)
end

function QuantifyPeaks(signals, f)
    numpeaks = 0
    signal = 0
    peaks = []
    indices = []
    newpeak = []
    for i in 1:length(signals)
        if signals[i] == 1
            if signal == 0
                numpeaks += 1
                push!(indices, i-1)
                if i > 1
                    newpeak = [0.0]
                end
            end
            push!(newpeak, f[i])
            signal = 1
        else
            if signal == 1
                push!(newpeak, 0)
                push!(peaks, newpeak)
            end
            signal = 0
        end
    end
    numpeaks, peaks, indices
end

optimise()
