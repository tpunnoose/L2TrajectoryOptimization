function CWDynamics!(ẋ, x, u, prob)
    # Linear dynamics model of the ego and target satellite.
    # Implements the Clohessy-Wiltshire equations
    # with control added as force along the (x,y,z) axis defined in
    # R. Wiltshire and W. Clohessy,
    # “Terminal guidance system for satellite rendezvous”.

    m_ego = prob.m_ego
    n_ = prob.n_
    ẋ[1:3] = x[4:6]
    ẋ[4] =  2*n_*x[5] + u[1]/m_ego
    ẋ[5] = -2*n_*x[4] + 3*n_^2*x[2] + u[2]/m_ego
    ẋ[6] = -n_^2*x[3] + u[3]/m_ego
end
̇̇̇
