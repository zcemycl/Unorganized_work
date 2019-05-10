function new_ref = motion(phiI, phiM, phiIM, phiMF, mode, L, lineObj, xref, yref, n)
if mode == 'r'
    for i = 1:n
        Phi = phiI + phiIM*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi, L, 'r');
        animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'r');
    end
    
    for i = 1:n
        Phi2 = phiM + phiMF*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi2, L, 'r');
        [new_ref] = animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'r');
    end
elseif mode == 'f'
    for i = 1:n
        Phi = phiI + phiIM*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi, L, 'f');
        animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'f');
    end
    
    for i = 1:n
        Phi2 = phiM + phiMF*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi2, L, 'f');
        [new_ref] = animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'f');
    end
end
end