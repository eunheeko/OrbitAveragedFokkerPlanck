## Local first-order diffusion coefficient in velocitty change parallel to v: (δv∥)
function cal_parl1(v::Float64, bubparam::BubbleParameters)
    G = bubparam.G
    ρBUB = bubparam.ρBUB_IN
    RBUB = bubparam.RBUB_IN
    TBUB = bubparam.TBUB
    dnBUB = bubparam.dnBUB
    
    # C = -10752 / 25 * sqrt(14/5) * pi^(3) * dnBUB * G^2 * ρBUB^2 * RBUB^8 * TBUB^3 # constant
    C = -43008 / 25 * sqrt(14/5) * pi^(4) * dnBUB * G^2 * ρBUB^2 * RBUB^8 * TBUB^3 # constant

    parl1 = C / (100 * pi * RBUB^2 + 7 * v^2 * TBUB^2)^(5/2) * v
    
    return parl1
end

## Local first-order diffusion coefficient in velocitty change perpendicular to v: (δv⟂)
function cal_perp1(v::Float64, bubparam::BubbleParameters)
    G = bubparam.G
    ρBUB = bubparam.ρBUB_IN
    RBUB = bubparam.RBUB_IN
    TBUB = bubparam.TBUB
    dnBUB = bubparam.dnBUB
    
    # C = -1792 / 125 * sqrt(14/5) * pi^(2) * dnBUB * G^2 * ρBUB^2 * RBUB^6 * TBUB^3 # constant
    C = -7168 / 125 * sqrt(14/5) * pi^(3) * dnBUB * G^2 * ρBUB^2 * RBUB^6 * TBUB^3 # constant

    perp1 = C / (100 * pi * RBUB^2 + 7 * v^2 * TBUB^2)^(3/2) * v
    
    return perp1
end

## Local second-order diffusion coefficient in velocitty change parallel to v: (δv∥)²
function cal_parl2(v::Float64, bubparam::BubbleParameters)
    G = bubparam.G
    ρBUB = bubparam.ρBUB_IN
    RBUB = bubparam.RBUB_IN
    TBUB = bubparam.TBUB
    dnBUB = bubparam.dnBUB
    
    # C = 512 / 25 * sqrt(14/5) * pi^(3) * dnBUB * G^2 * ρBUB^2 * RBUB^8 * TBUB^1 # constant
    C = 4096 / 25 * sqrt(14/5) * pi^(4) * dnBUB * G^2 * ρBUB^2 * RBUB^8 * TBUB^1 # constant

    parl2 = C / (100 * pi * RBUB^2 + 7 * v^2 * TBUB^2)^(3/2)
    
    return parl2
end

## Local second-order diffusion coefficient in velocity change perpendicular to v: (δv⟂)² 
function cal_perp2(v::Float64, bubparam::BubbleParameters)
    G = bubparam.G
    ρBUB = bubparam.ρBUB_IN
    RBUB = bubparam.RBUB_IN
    TBUB = bubparam.TBUB
    dnBUB = bubparam.dnBUB

    # C = 256 / 625 * sqrt(14/5) * pi^(2) * dnBUB * G^2 * ρBUB^2 * RBUB^6 * TBUB^1 # constant
    C = 2048 / 625 * sqrt(14/5) * pi^(3) * dnBUB * G^2 * ρBUB^2 * RBUB^6 * TBUB^1 # constant

    perp2 = C / (100 * pi * RBUB^2 + 7 * v^2 * TBUB^2)^(1/2) 
    
    return perp2
end

## Local first-order diffusion coefficient in energy: (δE)
function cal_dE(r::Float64, v::Float64, bubparam::BubbleParameters)

    parl1 = cal_parl2(v, bubparam)
    parl2 = cal_parl2(v, bubparam)
    perp2 = cal_perp2(v, bubparam)

    dE = v * parl1 + 0.5 * (parl2 + perp2)

    return dE
end


## Local first-order diffusion coefficient in angular momentum: (δL)
function cal_dL(r::Float64, v::Float64, L::Float64, bubparam::BubbleParameters)

    parl1 = cal_parl2(v, bubparam)
    perp2 = cal_perp2(v, bubparam)

    dL = (L / v) * parl1 + (r^2 / (4 * L)) * perp2

    return dL
end


## Local second-order diffusion coefficient in angular momentum: (δL)² 
function cal_dLL(r::Float64, v::Float64, L::Float64, bubparam::BubbleParameters)

    parl2 = cal_parl2(v, bubparam)
    perp2 = cal_perp2(v, bubparam)

    dLL = (L^2 / v^2) * parl2 + 0.5 * (r^2 - (L^2 / v^2)) * perp2

    return dLL
end

## Local second-order diffusion coefficient in energy: (δE)² 
function cal_dEE(r::Float64, v::Float64, bubparam::BubbleParameters)

    parl2 = cal_parl2(v, bubparam)

    dEE = v^2 * parl2

    return dEE
end

## Local second-order diffusion coefficient in angular momentum: (δEδL)
function cal_dEL(r::Float64, v::Float64, L::Float64, bubparam::BubbleParameters)

    parl2 = cal_parl2(v, bubparam)

    dEL = L * parl2

    return dEL
end