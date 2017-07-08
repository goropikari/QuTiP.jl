n_vec = [1, 2, 3]
pbar = TextProgressBar(3)
@test for n in n_vec
    pbar[:update](n)
    exp(n)
    pbar[:finished]()
end == nothing
