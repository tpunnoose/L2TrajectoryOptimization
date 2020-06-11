function SelectControl(k)
	return (3*(k-1)+1):(3*(k-1)+3)
end

function SelectState(k)
	return (6*(k-1)+1):(6*(k-1)+6)
end

function StateTrajectoryToArray(X, params)
	reshape(X, params.n, :)
end

function ControlTrajectoryToArray(U, params)
	reshape(U, params.m, :)
end
