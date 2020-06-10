function ADMM(parameters, num_iter)
	n = 6
	N = 10
	m = 3

	X = zeros(n*(N+1))
	U = zeros(m*N)
	Y = zeros(m*N)
	Λ̄ = zeros(m*N)

	for i=1:num_iter
		(X, U) = LQR(ρ, Y, Λ̄)
		for k=1:N
			β = α/ρ
			v = U[select3(k)] + Λ̄[select3(k)]
			Y[select3(k)] = l2_prox(β, v)

			Λ̄[select3(k)] += U[select3(k)] - Y[select3(k)]
		end
	end
end

function l2_prox(β, v)
	return max(0, 1-β/norm(v, 2))*v
end

function select3(k)
	return (3*(i-1)+1):(3*(i-1)+3)
end
