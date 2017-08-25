# pretty-printing
# Qobjs are displayed by mathjax on IJulia
Base.show(io::IO, ::MIME"text/latex", s::Quantum) =  begin
    if s[:__class__][:__name__] == "Qobj"
        print(io, s[:_repr_latex_]())
    else
        print(s.o)
    end
end

function Base.show(io::IO, ::MIME"text/latex", s::Vector{Quantum})
    for item in s
        show(io, "text/latex", item)
    end
end
