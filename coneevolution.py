def cone(D,age):
# outputs: 
# topo -- topographic elevation of the cone [m]
# distance -- distance along the cone profile [m]
# coneslope -- mean slope of cone [deg]
# inputs: 
# D -- landscape diffusivity [m^2/kyr]
# age -- age of the cone [kyr]

    import numpy as np

    # Topo and Grid Parameters
    conediam=876
    craterdiam=292
    nx=600
    dx=5
    coneslope=0.575

    # Create Initial Cone Topography
    icenter=round(nx/2)
    cone_radius=conediam/2
    crater_radius=craterdiam/2
    topo=1000+np.zeros(nx)
    dist=np.zeros(nx)
    distance=np.zeros(nx)
    for i in range(nx):
        distance[i]=i*dx
        dist[i]=abs(i*dx-icenter*dx);
        if dist[i]<=crater_radius:
            topo[i]=topo[i]+dist[i]*coneslope;
        if ((dist[i]>crater_radius) and (dist[i]<cone_radius)):
            topo[i]=1000+crater_radius*coneslope-(dist[i]-crater_radius)*coneslope;
        if dist[i]>cone_radius:
            topo[i]=1000+crater_radius*coneslope-(cone_radius-crater_radius)*coneslope;

    topo=topo-min(topo);

    # 1d Diffusion
    tend=age           # kyr

    toponew=np.zeros(nx)
    dt=0.1*dx*dx/D
    slope=np.zeros(nx)

    t=0
    while t<tend:
    
        for i in range(1,nx-1):
            toponew[i]=topo[i]+D*dt/(dx*dx)*(topo[i+1]-2*topo[i]+topo[i-1])

        toponew[0]=topo[0]+D*dt/(dx*dx)*(topo[1]-2*topo[0]+topo[0])
        toponew[nx-1]=topo[nx-1]+D*dt/(dx*dx)*(topo[nx-1]-2*topo[nx-1]+topo[nx-2])
    
        t+=dt
        topo=toponew[:]
        
    slope=np.gradient(topo,dx)
    slope=180/np.pi*np.arctan(abs(slope))

    maxtopo=max(topo)

    # Find array index associated with lowest elevation on the cone    
    coneedge=-1
    i=0
    while coneedge==-1:
        if slope[i]>3:
            coneedge=i
        i+=1
        
    # Find array index associated with highest elevation on the cone
    maxind=-1;
    i=0;
    while maxind==-1:
        if ((topo[i+1]<topo[i]) or (i==nx-1)):
            maxind=i
        i+=1

    coneslope=(maxtopo-topo[coneedge])/(dx*(maxind-coneedge))
    coneslope=180/np.pi*np.arctan(coneslope)

    return (topo,distance,coneslope);

