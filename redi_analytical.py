def redi_simple(lat,kappa,v,La,modelh,ntime,dt,Csouth,Cnorth):
    from scipy import signal
    import numpy as np
    ny=144
    C=np.ones((ny,40000))*3
    Cnew=np.ones(ny)*3
    Cold=np.ones(ny)*3
    Cs=np.zeros(ny)
    Cs[:]= 20*(((np.cos(np.pi*lat[:]/140))**2) + (0.1*np.exp(- (lat[:]/120-1)**2)))

    dy=108086
    #note that walt has units m/s because it is actually w
    watl=-np.diff(v)/dy/La[:-1]
    #v has units m^3/s

    kappav=2*10**-5

    timestep=dt*24*3600*365
    for time in range(1, ntime):
        Cnew[1:-1]=(Cold[1:-1]+
                timestep*((kappa[2:]*La[2:]*Cold[2:]*(modelh[1:-1]+modelh[2:])/2
                           -kappa[2:]*La[2:]*Cold[1:-1]*(modelh[1:-1]+modelh[2:])/2- 
                           kappa[1:-1]*La[1:-1]*Cold[1:-1]*(modelh[1:-1]+modelh[0:-2])/2
                           +kappa[1:-1]*La[1:-1]*Cold[0:-2]*(modelh[1:-1]+modelh[0:-2])/2)/dy**2/La[1:-1]/modelh[1:-1]
                          -(v[2:]*(Cold[2:]+Cold[1:-1])-v[1:-1]*(Cold[1:-1]+Cold[0:-2]))/dy/2/La[1:-1]/modelh[1:-1]
                           -watl[1:]*Cold[1:-1]/modelh[1:-1]+kappav*(Cs[1:-1]-Cold[1:-1])/modelh[1:-1]/modelh[1:-1]))#
        #Boundary conditions are set here. 
        Cnew[0]=Csouth
        Cnew[-1]=Cnorth
        Cold=Cnew
        timeend=time-ntime+40000
        if timeend>=0:
            C[:,timeend]=Cnew
    return C


