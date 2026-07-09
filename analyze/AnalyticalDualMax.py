import numpy as np

def multipole(A1, A2, A3, r0=0):
    def func(r, th, t):
        # Frequently occurring exponentials
        Emm = np.exp(-(r - r0 - t)**2)
        Epm = np.exp(-(r + r0 - t)**2)
        Emp = np.exp(-(r - r0 + t)**2)
        Epp = np.exp(-(r + r0 + t)**2)

        # Useful coordinate combinations
        xm = r - r0 - t
        xp = r + r0 - t
        ym = r - r0 + t
        yp = r + r0 + t

        # ------------------------------------------------------------------
        # Common combinations
        # ------------------------------------------------------------------

        S0 = 0.5*(Emm + Epm)*(r - t) + 0.5*(Emp + Epp)*(r + t)

        S1 = (
            0.5*(Emm + Epm)
            + 0.5*(Emp + Epp)
            + 0.5*(-2*Emm*xm - 2*Epm*xp)*(r - t)
            + 0.5*(r + t)*(-2*Emp*ym - 2*Epp*yp)
        )

        S2 = (
            0.5*(-2*Emm - 2*Epm + 4*Emm*xm**2 + 4*Epm*xp**2)*(r - t)
            - 2*Emm*xm
            - 2*Epm*xp
            - 2*Emp*ym
            - 2*Epp*yp
            + 0.5*(r + t)*(
                -2*Emp
                -2*Epp
                +4*Emp*ym**2
                +4*Epp*yp**2
            )
        )

        S3 = (
            -2*Emm
            -2*Epm
            -2*Emp
            -2*Epp

            +0.5*(-2*Emm -2*Epm +4*Emm*xm**2 +4*Epm*xp**2)

            +0.5*(
                12*Emm*xm
                -8*Emm*xm**3
                +12*Epm*xp
                -8*Epm*xp**3
            )*(r - t)

            +4*Emm*xm**2
            +4*Epm*xp**2
            +4*Emp*ym**2
            +4*Epp*yp**2

            +0.5*(-2*Emp -2*Epp +4*Emp*ym**2 +4*Epp*yp**2)

            +0.5*(r + t)*(
                12*Emp*ym
                -8*Emp*ym**3
                +12*Epp*yp
                -8*Epp*yp**3
            )
        )

        # ------------------------------------------------------------------
        # Final expression
        # ------------------------------------------------------------------

        expr = (
            0.5*A1*(
                S0/r**2
                - S1/r
            )*np.sin(th)

            +0.5*A2*(
                S0/r**3
                - S1/r**2
                + S2/(3*r)
            )*np.cos(th)*np.sin(th)

            +0.5*A3*(
                S0/r**4
                - S1/r**3
                + (2*S2)/(5*r**2)
                - S3/(15*r)
            )*(5*np.cos(th)**2 - 1)*np.sin(th)
        )
        return expr
    return func