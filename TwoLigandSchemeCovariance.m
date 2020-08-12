\[CapitalSigma]data = 
    {{(kp*((c1*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)^2)/
          (kf*(1 + koff/kf)) + (c2*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)^
            2)/(kf2*(1 + koff2/kf2)) + 
         (2*((c1*c2^2*koff2*(1 + koff2/kf2)*kon^3)/(kf2^2*koff) + 
            (c2*koff^4*kon*(koff2^3/kf2^3 + (c1^2*kon^2)/koff^2 + 
               (c1*koff2*kon*(-1 + (2*c1*kon)/koff + (c2*kon)/koff2))/
                (kf2*koff) + (koff2^2*(1 + (c2*kon)/koff2 + (c1^2*kon^2)/
                   koff^2 + (c2^2*kon^2)/koff2^2 + (c1*kon*(-1 + (2*c2*kon)/
                      koff2))/koff))/kf2^2))/(kf^3*koff2) + 
            (c2*koff*koff2*kon*(kf*(1 + (c2*kon)/koff2 + (c1^2*kon^2)/
                  koff^2 + (c2^2*kon^2)/koff2^2 + (koff2*(1 + (c1*kon)/koff)^
                    2)/kf2 + (c1*kon*(2 + (c2*kon)/koff2))/koff) + 
               (c1*kf2*(1 + koff2/kf2)*kon*(-1 + (c1*kon)/koff - (c2*kon)/
                   koff2 + (koff2*(-1 + (c1*kon)/koff + (2*c2*kon)/koff2))/
                   kf2))/koff))/(kf*kf2^2) + (koff^3*((c1*koff2^4*kon)/
                (kf2^3*koff) + (c1^2*c2*kf*kon^3)/(koff^2*koff2) + 
               (koff2^3*((c2*kf*kon*(3 + (2*c1*kon)/koff))/koff2 + 
                  (c1*kf2*kon*(3 + (2*c2*kon)/koff2))/koff))/kf2^3 + 
               (c1*koff2*kon*(kf2*(1 + (c2*kon)/koff2)^2 + (c2*kf*kon*
                    (-2 + (c1*kon)/koff + (2*c2*kon)/koff2))/koff2))/
                (kf2*koff) + (koff2^2*((c1*kf2*kon*(3 + (4*c2*kon)/koff2 + 
                     (c2^2*kon^2)/koff2^2))/koff + (c2*kf*kon*(3 + (3*c2*kon)/
                      koff2 + (3*c2^2*kon^2)/koff2^2 + (5*c1*c2*kon^2)/
                      (koff*koff2)))/koff2))/kf2^2))/kf^3 + 
            (koff^2*koff2*((c1*kf2*(1 + koff2/kf2)*kon*((c1^2*kon^2)/koff^2 + 
                  (c1*kon*(1 + (c2*kon)/koff2))/koff + (1 + (c2*kon)/koff2)^
                   2 + (koff2^2*(1 + (c1*kon)/koff - (c2*kon)/koff2 + 
                     (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/koff2^2 + 
                     (2*c1*c2*kon^2)/(koff*koff2)))/kf2^2 + 
                  (koff2*(2 + (c2*kon)/koff2 + (2*c1^2*kon^2)/koff^2 - 
                     (c2^2*kon^2)/koff2^2 + (c1*kon*(2 + (3*c2*kon)/koff2))/
                      koff))/kf2))/koff + (c2*kf*kon*
                 ((c1*kon*(-1 - (c1*kon)/koff + (c2*kon)/koff2))/koff + 
                  (koff2^2*(3 + (4*c1*kon)/koff + (c1^2*kon^2)/koff^2))/
                   kf2^2 + (koff2*((c1*kon*(3 + (4*c2*kon)/koff2))/koff + 
                     3*(1 + (c2*kon)/koff2 + (c2^2*kon^2)/koff2^2)))/kf2))/
                koff2))/(kf^2*kf2))*kp)/(koff*(1 + koff/kf)^3*koff2*
           (1 + koff2/kf2)^3))*t)/(1 + (c1*kon)/koff + (c2*kon)/koff2)^3, 
      -(kp*t*((2*(1 + koff/kf)*(1 + koff2/kf2)*(c1*kon + c2*kon)*
            (-((c1*c2*koff2*kon^2)/(kf2*koff)) + (c2*koff^2*kon*(koff2/kf2 - 
                (c1*kon)/koff))/(kf*koff2) + (koff*koff2*(
                (c2*kf*kon*(1 + (c1*kon)/koff))/koff2 + (c1*kf2*kon*
                  (1 + koff2/kf2 + (c2*kon)/koff2))/koff))/(kf*kf2)))/
           (koff*koff2) + (8*(1 + (c1*kon)/koff + (c2*kon)/koff2)*
            (-koff2^2 - 2*c2*koff2*kon - c2^2*kon^2 - 
             koff^2*(1 + (c1*kon)/koff)^2 - 2*kf*Sqrt[koff^2*
                 (1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 
                2*koff*koff2*(-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/
                      koff2))/koff)] + koff2*Sqrt[koff^2*(1 + (c1*kon)/koff)^
                  2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                 (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                   koff)] + c2*kon*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
                koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - 
                  (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)] + 
             koff*(2*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2 - 
                 (c1*c2*kon^2)/(koff*koff2)) + (-1 + (c1*kon)/koff)*
                Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                   (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                     koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)]))*
            ((c1*c2^2*koff2*(1 + koff2/kf2)*kon^3)/(kf2*koff) + 
             (c1*c2*koff^3*kon^2*(-(koff2^2/kf2^2) + (c1*kon)/koff + 
                (koff2*((2*c1*kon)/koff + (c2*kon)/koff2))/kf2))/
              (kf^2*koff2) + (c2*koff*koff2*kon*((c1*kf2*(1 + koff2/kf2)*kon*
                  ((2*c2*kon)/kf2 + (c1*kon)/koff + (c1*koff2*kon)/
                    (kf2*koff)))/koff + kf*(1 + (c1*kon)/koff + (c2*kon)/
                   koff2 + (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/koff2^2 + 
                  (koff2*(1 + (c1*kon)/koff + (c1^2*kon^2)/koff^2 - 
                     (c1*c2*kon^2)/(koff*koff2)))/kf2)))/(kf*kf2) + 
             (koff^3*((c1^2*c2*kf*kon^3)/(koff^2*koff2) + 
                (koff2^3*((c1*kon)/koff + (c2*kon)/koff2 + (2*c1*c2*kon^2)/
                    (koff*koff2)))/kf2^2 + (c1*koff2*kon*((2*c2*kf*kon*
                     ((c1*kon)/koff + (c2*kon)/koff2))/koff2 + 
                   kf2*(1 + (c2*kon)/koff2 + (c2^2*kon^2)/koff2^2 - 
                     (c1*c2*kon^2)/(koff*koff2))))/(kf2*koff) + 
                (koff2^2*(-((c1*c2*kf*kon^2*(2 + (c1*kon)/koff))/(koff*
                      koff2)) + kf2*((2*c1*kon)/koff + (c2*kon)/koff2 + 
                     (c2^2*kon^2)/koff2^2 + (3*c1*c2*kon^2)/(koff*koff2) + 
                     (c2^3*kon^3)/koff2^3 + (2*c1*c2^2*kon^3)/(koff*
                       koff2^2) - (c1^2*c2*kon^3)/(koff^2*koff2))))/kf2^2))/
              kf^2 + (koff^2*koff2*(-((c1*c2*kf2*(1 + koff2/kf2)*kon^2*
                   (1 + koff2/kf2 + (c2*kon)/koff2))/koff) - 
                (c1*c2*kf^2*kon^2*(koff2/kf2 - (c2*kon)/koff2 + 
                   (c1*koff2*kon)/(kf2*koff)))/(koff*koff2) + kf*kf2*
                 ((c1^3*(1 + koff2/kf2)^2*kon^3)/koff^3 + 
                  (c1^2*(1 + koff2/kf2)*kon^2*(1 + koff2/kf2 + (2*c2*kon)/
                      kf2))/koff^2 + (2*c2*kon*(1 + koff2/kf2 + (c2*kon)/
                      koff2 + (c2^2*kon^2)/koff2^2))/kf2 + 
                  (c1*kon*(1 + (c2*kon)/koff2 + (c2^2*kon^2)/koff2^2 + 
                     (2*koff2*(1 + (c2*kon)/koff2)^2)/kf2 + (koff2^2*
                       (1 + (3*c2*kon)/koff2 - (c2^2*kon^2)/koff2^2))/kf2^2))/
                   koff)))/(kf^2*kf2)))/(Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
              koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - 
                (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)]*
            (koff2 + c2*kon + koff*(1 + (c1*kon)/koff) - 
             Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                 2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                 (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*(-koff2 - c2*kon + 
             kf*(2 + koff/kf - (c1*kon)/kf) + Sqrt[koff^2*(1 + (c1*kon)/koff)^
                 2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*
            (koff2 + c2*kon + koff*(1 + (c1*kon)/koff) + 
             Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                 2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                 (c1*kon*(-1 + (c2*kon)/koff2))/koff)])) + 
          (1 + koff/kf)*(1 + koff2/kf2)*(c1*kon + c2*kon)*
           (1 + (c1*kon)/koff + (c2*kon)/koff2)*((c2*kon)/kf2 + 
            (koff*((c2*kon)/kf2 + (c1*kon)/koff + (c1*koff2*kon)/(kf2*koff)))/
             kf)*t - (4*koff*(1 + koff/kf)*koff2*(1 + koff2/kf2)*
            (c1*kon + c2*kon)*(1 + (c1*kon)/koff + (c2*kon)/koff2)^2*
            ((c2*kon)/kf2 + (koff*((c2*kon)/kf2 + (c1*kon)/koff + 
                (c1*koff2*kon)/(kf2*koff)))/kf)*(koff2^2 + 2*c2*koff2*kon + 
             c2^2*kon^2 + koff^2*(1 + (c1*kon)/koff)^2 + 
             2*kf*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                 (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                  (c1*kon*(-1 + (c2*kon)/koff2))/koff)] - 
             koff2*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                 (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                  (c1*kon*(-1 + (c2*kon)/koff2))/koff)] - 
             c2*kon*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                 (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                  (c1*kon*(-1 + (c2*kon)/koff2))/koff)] - 
             koff*(2*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2 - 
                 (c1*c2*kon^2)/(koff*koff2)) + (-1 + (c1*kon)/koff)*
                Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                   (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                     koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)]))*t)/
           (Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                (c1*kon*(-1 + (c2*kon)/koff2))/koff)]*(koff2 + c2*kon + 
             koff*(1 + (c1*kon)/koff) - Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
               koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*
            (-koff2 - c2*kon + kf*(2 + koff/kf - (c1*kon)/kf) + 
             Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                 2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                 (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*(koff2 + c2*kon + 
             koff*(1 + (c1*kon)/koff) + Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
               koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                  koff)]))))/(2*(1 + koff/kf)^2*(1 + koff2/kf2)^2*
        (1 + (c1*kon)/koff + (c2*kon)/koff2)^3), 
      -((((2*c1*c2*koff2*(1 + koff2/kf2)*kon^2*(1 + koff2/kf2 + 
             (c2*kon)/kf2))/(kf2*koff) + (c2*koff^4*kon*
            ((2*c1*kon*(1 + (c1*kon)/koff))/koff + 
             (koff2^2*(-1 + (c1*kon)/koff + (3*c2*kon)/koff2 + (2*c1^2*kon^2)/
                 koff^2 + (2*c2^2*kon^2)/koff2^2 + (4*c1*c2*kon^2)/
                 (koff*koff2)))/kf2^2 + (koff2*(-1 + (c2*kon)/koff2 + 
                (4*c1^2*kon^2)/koff^2 + (c1*kon*(3 + (2*c2*kon)/koff2))/
                 koff))/kf2))/(kf^3*koff2) + 
          (koff*koff2*((c2*kf*kon*(-1 - (c1*kon)/koff + (c2*kon)/koff2 + 
                (2*c1*koff2^2*kon*(1 + (c1*kon)/koff))/(kf2^2*koff) + 
                (koff2*(-1 + (c1*kon)/koff + (3*c2*kon)/koff2 + 
                   (2*c1^2*kon^2)/koff^2 + (2*c2^2*kon^2)/koff2^2 + 
                   (2*c1*c2*kon^2)/(koff*koff2)))/kf2))/koff2 + 
             (c1*kf2*(1 + koff2/kf2)*kon*(-1 + (c1*kon)/koff - (c2*kon)/
                 koff2 + (2*koff2*(-1 + (c1*kon)/koff + (c2*kon)/koff2 - 
                   (c2^2*kon^2)/koff2^2 + (c1*c2*kon^2)/(koff*koff2)))/kf2 + 
                (koff2^2*(-1 + (c1*kon)/koff + (3*c2*kon)/koff2 + 
                   (4*c2^2*kon^2)/koff2^2 + (2*c1*c2*kon^2)/(koff*koff2)))/
                 kf2^2))/koff))/(kf*kf2) + 
          (koff^2*((2*c1*c2*kf*kon^2)/(koff*koff2) + (c1*koff2^4*kon*(-1 + 
                (c2*kon)/koff2 + (2*c1^2*kon^2)/koff^2 + (2*c2^2*kon^2)/
                 koff2^2 + (c1*kon*(3 + (4*c2*kon)/koff2))/koff))/
              (kf2^3*koff) + (koff2*((c2*kf*kon*(-3 + (c1*kon)/koff + 
                   (3*c2*kon)/koff2 - (2*c1^2*kon^2)/koff^2 + (2*c1*c2*kon^2)/
                    (koff*koff2)))/koff2 + (c1*kf2*kon*(-1 + (3*c1*kon)/
                    koff + (c2*kon)/koff2 + (2*c1^2*kon^2)/koff^2 + 
                   (2*c2^2*kon^2)/koff2^2 + (2*c1*c2*kon^2)/(koff*koff2)))/
                 koff))/kf2 + (koff2^2*((c1*kf2*kon*(-3 + (9*c1*kon)/koff + 
                   (3*c2*kon)/koff2 + (6*c1^2*kon^2)/koff^2 + (8*c1*c2*kon^2)/
                    (koff*koff2)))/koff + (c2*kf*kon*(-3 + (3*c1*kon)/koff + 
                   (9*c2*kon)/koff2 + (6*c2^2*kon^2)/koff2^2 + 
                   (8*c1*c2*kon^2)/(koff*koff2)))/koff2))/kf2^2 + 
             (c1*koff2^3*kon*((2*c2*kf*kon*(2 + (c1*kon)/koff))/koff2 + 
                kf2*(-3 + (9*c1*kon)/koff + (3*c2*kon)/koff2 + (6*c1^2*kon^2)/
                   koff^2 + (10*c1*c2*kon^2)/(koff*koff2))))/(kf2^3*koff)))/
           kf^2 + (c2*koff^3*kon*((2*c1*(kf + kf2)*koff2^3*kon)/
              (kf2^3*koff) + (2*c1*kf*kon*(2 + (c1*kon)/koff))/koff + 
             (koff2*((2*c1*kf2*kon*(1 + (c2*kon)/koff2))/koff + 
                kf*(-3 + (5*c1*kon)/koff + (3*c2*kon)/koff2 + (2*c1^2*kon^2)/
                   koff^2 + (4*c1*c2*kon^2)/(koff*koff2))))/kf2 + 
             (koff2^2*((2*c1*kf2*kon*(2 + (c2*kon)/koff2))/koff + 
                kf*(-3 + (3*c1*kon)/koff + (9*c2*kon)/koff2 + (6*c2^2*kon^2)/
                   koff2^2 + (10*c1*c2*kon^2)/(koff*koff2))))/kf2^2))/
           (kf^3*koff2))*kp^2*t)/(koff*(1 + koff/kf)^3*koff2*
         (1 + koff2/kf2)^3*(1 + (c1*kon)/koff + (c2*kon)/koff2)^3)), 
      (kp*t*((2*c2*(1 + koff/kf)^2*kon*((c1*c2*koff2*(1 + koff2/kf2)*kon^2)/
             (kf2*koff) + (koff^2*((c1*c2*kon^2)/(koff*koff2) + (koff2^2*
                 (1 + (c1*kon)/koff))/kf2^2 + (koff2*(1 + (c1*kon)/koff + 
                  (c2*kon)/koff2 + (c2^2*kon^2)/koff2^2 + (2*c1*c2*kon^2)/
                   (koff*koff2)))/kf2))/kf + 
            (koff*koff2*(-((c1*kf2*(1 + koff2/kf2)*kon*(1 + koff2/kf2 + 
                   (c2*kon)/koff2))/koff) + kf*(1 + (c2*kon)/koff2 + 
                 (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/koff2^2 + 
                 (koff2*(1 + (c1*kon)/koff)^2)/kf2 + (c1*kon*(2 + (c2*kon)/
                     koff2))/koff)))/(kf*kf2)))/(koff*koff2) + 
         (2*c1*(1 + koff2/kf2)^2*kon*((c1*c2*koff2*kon^2)/(kf2*koff) + 
            (c2*koff^3*kon*(-(koff2/kf2) + (c1*kon)/koff))/(kf^2*koff2) + 
            (koff^2*((c1*c2*kf*kon^2)/(koff*koff2) + (koff2^2*
                 (1 + (c2*kon)/koff2))/kf2 + (koff2*(-((c2*kf*kon*
                     (2 + (c1*kon)/koff))/koff2) + kf2*(1 + (c2*kon)/koff2)^
                    2))/kf2))/kf^2 + (koff*koff2*(-((c2*kf*kon*
                  (1 + (c1*kon)/koff))/koff2) + kf2*((c1^2*kon^2)/koff^2 + 
                 (c1*kon*(1 + (c2*kon)/koff2))/koff + (1 + (c2*kon)/koff2)^
                  2 + (koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2 + 
                    (c1^2*kon^2)/koff^2 + (2*c1*c2*kon^2)/(koff*koff2)))/
                  kf2)))/(kf*kf2)))/(koff*koff2) - 
         (2*c1*(1 + koff2/kf2)^2*kon*(-((c1*c2*koff2*kon^2)/koff) + 
            (c2*koff^3*kon*(-((c1*kf*kon)/koff) + (koff2*(kf2 - (c1*kf*kon)/
                   koff + (c1*kf2*kon)/koff))/kf2))/(kf^2*koff2) + 
            (koff*koff2*(c2*kon*(1 - (c1*kon)/koff + (c2*kon)/koff2) + kf*
                ((c2*kon)/koff2 + (c1*kon*(1 + koff2/kf2 + (c2*kon)/kf2 + 
                    (2*c2*kon)/koff2))/koff)))/kf + 
            (koff^2*(-((c1*c2*kf^2*kon^2)/(koff*koff2)) + (koff2^2*
                 ((c2*kf2*kon*(1 + (c2*kon)/koff2))/koff2 + (c1*kf*kon*
                    (2 + (c1*kon)/koff + (2*c2*kon)/koff2))/koff))/kf2 + 
               (kf*koff2*(-((c1*c2*kf*kon^2)/(koff*koff2)) + 
                  kf2*((2*c1*kon)/koff + (2*c2*kon)/koff2 + (c1^2*kon^2)/
                     koff^2 + (4*c1*c2*kon^2)/(koff*koff2))))/kf2))/kf^2))/
          (kf*koff*koff2) - (2*c2*(1 + koff/kf)^2*kon*
           (-((c1*c2*koff2*(1 + koff2/kf2)*kon^2)/koff) + 
            (koff^2*(-((c1*c2*kf*kon^2)/(koff*koff2)) + (koff2*
                 ((c2*kf2*kon*(1 + (c1*kon)/koff))/koff2 + (c1*kf*kon*
                    (1 + (c1*kon)/koff - (c2*kon)/koff2))/koff))/kf2 + 
               (koff2^2*((c1*kf*kon*(1 + (c1*kon)/koff))/koff + 
                  (c2*kf2*kon*(2 + (2*c1*kon)/koff + (c2*kon)/koff2))/koff2))/
                kf2^2))/kf + (koff*koff2*(-((c1*c2*(1 + koff2/kf2)*kon^2)/
                 koff) + kf*((c2*kon*(1 + (koff2*(2 + (c2*kon)/koff2))/kf2))/
                  koff2 + (c1*kon*(1 + (2*c2*kon)/koff2 + (koff2^2*
                      (1 + (c2*kon)/koff2))/kf2^2 + (koff2*(2 + (4*c2*kon)/
                        koff2))/kf2))/koff)))/kf))/(kf2*koff*koff2) + 
         c2*(1 + koff/kf)^2*(1 + koff2/kf2)*kon*(1 + (c1*kon)/koff + 
           (c2*kon)/koff2)*((c2*kon)/kf2 + 
           (koff*((c2*kon)/kf2 + (c1*kon)/koff + (c1*koff2*kon)/(kf2*koff)))/
            kf)*t + c1*(1 + koff/kf)*(1 + koff2/kf2)^2*kon*
          (1 + (c1*kon)/koff + (c2*kon)/koff2)*((c2*kon)/kf2 + 
           (koff*((c2*kon)/kf2 + (c1*kon)/koff + (c1*koff2*kon)/(kf2*koff)))/
            kf)*t + (c2*(1 + koff/kf)^2*(1 + koff2/kf2)*kon*
           (1 + (c1*kon)/koff + (c2*kon)/koff2)*(c2*kon + 
            (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*kon)/koff))/kf)*t)/kf2 + 
         (c1*(1 + koff/kf)*(1 + koff2/kf2)^2*kon*(1 + (c1*kon)/koff + 
            (c2*kon)/koff2)*(c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*
                 kon)/koff))/kf)*t)/kf - 2*(1 + koff/kf)*(1 + koff2/kf2)*
          (1 + (c1*kon)/koff + (c2*kon)/koff2)*((c2*kon)/kf2 + 
           (koff*((c2*kon)/kf2 + (c1*kon)/koff + (c1*koff2*kon)/(kf2*koff)))/
            kf)*(c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*kon)/koff))/
            kf)*t))/(2*(1 + koff/kf)^3*(1 + koff2/kf2)^3*
        (1 + (c1*kon)/koff + (c2*kon)/koff2)^3)}, 
     {-(kp*t*((2*(1 + koff/kf)*(1 + koff2/kf2)*(c1*kon + c2*kon)*
            (-((c1*c2*koff2*kon^2)/(kf2*koff)) + (c2*koff^2*kon*(koff2/kf2 - 
                (c1*kon)/koff))/(kf*koff2) + (koff*koff2*(
                (c2*kf*kon*(1 + (c1*kon)/koff))/koff2 + (c1*kf2*kon*
                  (1 + koff2/kf2 + (c2*kon)/koff2))/koff))/(kf*kf2)))/
           (koff*koff2) + (8*(1 + (c1*kon)/koff + (c2*kon)/koff2)*
            (-koff2^2 - 2*c2*koff2*kon - c2^2*kon^2 - 
             koff^2*(1 + (c1*kon)/koff)^2 - 2*kf*Sqrt[koff^2*
                 (1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 
                2*koff*koff2*(-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/
                      koff2))/koff)] + koff2*Sqrt[koff^2*(1 + (c1*kon)/koff)^
                  2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                 (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                   koff)] + c2*kon*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
                koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - 
                  (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)] + 
             koff*(2*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2 - 
                 (c1*c2*kon^2)/(koff*koff2)) + (-1 + (c1*kon)/koff)*
                Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                   (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                     koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)]))*
            ((c1*c2^2*koff2*(1 + koff2/kf2)*kon^3)/(kf2*koff) + 
             (c1*c2*koff^3*kon^2*(-(koff2^2/kf2^2) + (c1*kon)/koff + 
                (koff2*((2*c1*kon)/koff + (c2*kon)/koff2))/kf2))/
              (kf^2*koff2) + (c2*koff*koff2*kon*((c1*kf2*(1 + koff2/kf2)*kon*
                  ((2*c2*kon)/kf2 + (c1*kon)/koff + (c1*koff2*kon)/
                    (kf2*koff)))/koff + kf*(1 + (c1*kon)/koff + (c2*kon)/
                   koff2 + (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/koff2^2 + 
                  (koff2*(1 + (c1*kon)/koff + (c1^2*kon^2)/koff^2 - 
                     (c1*c2*kon^2)/(koff*koff2)))/kf2)))/(kf*kf2) + 
             (koff^3*((c1^2*c2*kf*kon^3)/(koff^2*koff2) + 
                (koff2^3*((c1*kon)/koff + (c2*kon)/koff2 + (2*c1*c2*kon^2)/
                    (koff*koff2)))/kf2^2 + (c1*koff2*kon*((2*c2*kf*kon*
                     ((c1*kon)/koff + (c2*kon)/koff2))/koff2 + 
                   kf2*(1 + (c2*kon)/koff2 + (c2^2*kon^2)/koff2^2 - 
                     (c1*c2*kon^2)/(koff*koff2))))/(kf2*koff) + 
                (koff2^2*(-((c1*c2*kf*kon^2*(2 + (c1*kon)/koff))/(koff*
                      koff2)) + kf2*((2*c1*kon)/koff + (c2*kon)/koff2 + 
                     (c2^2*kon^2)/koff2^2 + (3*c1*c2*kon^2)/(koff*koff2) + 
                     (c2^3*kon^3)/koff2^3 + (2*c1*c2^2*kon^3)/(koff*
                       koff2^2) - (c1^2*c2*kon^3)/(koff^2*koff2))))/kf2^2))/
              kf^2 + (koff^2*koff2*(-((c1*c2*kf2*(1 + koff2/kf2)*kon^2*
                   (1 + koff2/kf2 + (c2*kon)/koff2))/koff) - 
                (c1*c2*kf^2*kon^2*(koff2/kf2 - (c2*kon)/koff2 + 
                   (c1*koff2*kon)/(kf2*koff)))/(koff*koff2) + kf*kf2*
                 ((c1^3*(1 + koff2/kf2)^2*kon^3)/koff^3 + 
                  (c1^2*(1 + koff2/kf2)*kon^2*(1 + koff2/kf2 + (2*c2*kon)/
                      kf2))/koff^2 + (2*c2*kon*(1 + koff2/kf2 + (c2*kon)/
                      koff2 + (c2^2*kon^2)/koff2^2))/kf2 + 
                  (c1*kon*(1 + (c2*kon)/koff2 + (c2^2*kon^2)/koff2^2 + 
                     (2*koff2*(1 + (c2*kon)/koff2)^2)/kf2 + (koff2^2*
                       (1 + (3*c2*kon)/koff2 - (c2^2*kon^2)/koff2^2))/kf2^2))/
                   koff)))/(kf^2*kf2)))/(Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
              koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - 
                (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)]*
            (koff2 + c2*kon + koff*(1 + (c1*kon)/koff) - 
             Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                 2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                 (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*(-koff2 - c2*kon + 
             kf*(2 + koff/kf - (c1*kon)/kf) + Sqrt[koff^2*(1 + (c1*kon)/koff)^
                 2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*
            (koff2 + c2*kon + koff*(1 + (c1*kon)/koff) + 
             Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                 2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                 (c1*kon*(-1 + (c2*kon)/koff2))/koff)])) + 
          (1 + koff/kf)*(1 + koff2/kf2)*(c1*kon + c2*kon)*
           (1 + (c1*kon)/koff + (c2*kon)/koff2)*((c2*kon)/kf2 + 
            (koff*((c2*kon)/kf2 + (c1*kon)/koff + (c1*koff2*kon)/(kf2*koff)))/
             kf)*t - (4*koff*(1 + koff/kf)*koff2*(1 + koff2/kf2)*
            (c1*kon + c2*kon)*(1 + (c1*kon)/koff + (c2*kon)/koff2)^2*
            ((c2*kon)/kf2 + (koff*((c2*kon)/kf2 + (c1*kon)/koff + 
                (c1*koff2*kon)/(kf2*koff)))/kf)*(koff2^2 + 2*c2*koff2*kon + 
             c2^2*kon^2 + koff^2*(1 + (c1*kon)/koff)^2 + 
             2*kf*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                 (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                  (c1*kon*(-1 + (c2*kon)/koff2))/koff)] - 
             koff2*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                 (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                  (c1*kon*(-1 + (c2*kon)/koff2))/koff)] - 
             c2*kon*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                 (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                  (c1*kon*(-1 + (c2*kon)/koff2))/koff)] - 
             koff*(2*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2 - 
                 (c1*c2*kon^2)/(koff*koff2)) + (-1 + (c1*kon)/koff)*
                Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                   (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                     koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)]))*t)/
           (Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                (c1*kon*(-1 + (c2*kon)/koff2))/koff)]*(koff2 + c2*kon + 
             koff*(1 + (c1*kon)/koff) - Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
               koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*
            (-koff2 - c2*kon + kf*(2 + koff/kf - (c1*kon)/kf) + 
             Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                 2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                 (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*(koff2 + c2*kon + 
             koff*(1 + (c1*kon)/koff) + Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
               koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                  koff)]))))/(2*(1 + koff/kf)^2*(1 + koff2/kf2)^2*
        (1 + (c1*kon)/koff + (c2*kon)/koff2)^3), 
      ((c1*kon + c2*kon)*((2*c1*c2*koff*kon^2)/koff2 + (2*c1*c2*koff2*kon^2)/
          koff + koff*koff2*(1 + (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/koff2^2 - 
           (2*c1*c2*kon^2)/(koff*koff2)))*t)/(koff*koff2*
        (1 + (c1*kon)/koff + (c2*kon)/koff2)^3), 
      -(((-2*c2*kon)/(kf2*(1 + koff2/kf2)^2) - (2*c2*kon)/
           (koff2*(1 + koff2/kf2)^2) - (2*c1*c2*kon^2)/
           (kf2*koff*(1 + koff2/kf2)^2) + (2*c1*c2*kon^2)/
           (koff2^2*(1 + koff2/kf2)^2) + (4*c1*c2*kon^2)/
           (kf2*koff2*(1 + koff2/kf2)^2) + (2*c2^2*kon^2)/
           (kf2*koff2*(1 + koff2/kf2)^2) - (2*c1*c2*kon^2)/
           (koff*koff2*(1 + koff2/kf2)^2) - (2*c1^2*c2*kon^3)/
           (kf2*koff^2*(1 + koff2/kf2)^2) + (2*c1*c2^2*kon^3)/
           (kf2*koff2^2*(1 + koff2/kf2)^2) + (2*c2^3*kon^3)/
           (kf2*koff2^2*(1 + koff2/kf2)^2) + (2*c1^2*c2*kon^3)/
           (koff*koff2^2*(1 + koff2/kf2)^2) + (2*c1*c2^2*kon^3)/
           (koff*koff2^2*(1 + koff2/kf2)^2) - (2*c1^2*c2*kon^3)/
           (koff^2*koff2*(1 + koff2/kf2)^2) + (4*c1^2*c2*kon^3)/
           (kf2*koff*koff2*(1 + koff2/kf2)^2) + (4*c1*c2^2*kon^3)/
           (kf2*koff*koff2*(1 + koff2/kf2)^2) - (2*c1*c2^2*kon^3)/
           (koff*koff2*(koff + (koff*koff2)/kf2)) + 
          (2*c1*(kf + koff - koff2)*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)^
             2*(c2*kon + kf*(-1 + (koff*(-1 + (c1*kon)/koff))/kf)))/
           (kf^2*(1 + koff/kf)^2*(kf*(1 + koff/kf)*(-1 + (c1*kon)/kf) + 
             koff2*(1 - (c1*kon)/kf + (c2*kon)/koff2))) + 
          (2*(c1*kon + c2*kon)*((c2*koff*kon*(1 + koff/kf + (c1*kon)/kf))/
              koff2 + (c1*koff2^2*kon*(1 + (c2*kon)/koff2))/(kf2*koff) - 
             (c1*koff2*kon*((c2*koff*kon)/koff2 + kf2*(-1 + (c2*koff*kon)/
                   (kf*koff2))))/(kf2*koff)))/(koff*(1 + koff/kf)*koff2*
            (1 + koff2/kf2)) + (2*kf*(1 + (c1*kon)/koff + (c2*kon)/koff2)*
            (koff^2*(1 + (c1*kon)/koff)^2 + kf2*(1 + koff2/kf2)*
              (koff2*(1 + (c2*kon)/koff2) + Sqrt[koff^2*(1 + (c1*kon)/koff)^
                   2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                  (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                    koff)]) - koff*(kf2*(1 + (c1*kon)/koff + 
                 (koff2*(2 + (c2*kon)/koff2 - (c1*kon*(-2 + (c2*kon)/koff2))/
                     koff))/kf2) + (1 + (c1*kon)/koff)*Sqrt[
                 koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                    2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + (c1*kon*
                      (-1 + (c2*kon)/koff2))/koff)]))*
            (-(c1*koff^2*kon*(1 + (c1*kon)/koff)^2) + 
             koff^2*(koff2*(1 + (4*c1*kon)/koff + (c1^3*kon^3)/koff^3 - 
                 (2*c1^2*kon^2*(-2 + (c2*kon)/koff2))/koff^2) - (c1*kon*
                 (1 + (c1*kon)/koff)*(-kf2 + Sqrt[koff^2*(1 + (c1*kon)/koff)^
                      2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                     (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                       koff)]))/koff) - kf2*koff*((koff2^2*(2 + (c2*kon)/
                   koff2 - (2*c1^2*kon^2*(-1 + (c2*kon)/koff2))/koff^2 + 
                  (c1*kon*(4 + (c2^2*kon^2)/koff2^2))/koff))/kf2 - 
               (c1*kon*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                    (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                      koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])/koff + 
               (koff2*(kf2*(1 + (c1*kon)/koff + (c1^2*kon^2)/koff^2 - 
                    (c1*c2*kon^2)/(koff*koff2)) - (1 + (c1^2*kon^2)/koff^2 - 
                    (c1*kon*(-2 + (c2*kon)/koff2))/koff)*Sqrt[
                    koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/
                         koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                       (c1*kon*(-1 + (c2*kon)/koff2))/koff)]))/kf2) + 
             kf2*koff2*((koff2^2*(1 + (c1*kon)/koff + (c2*kon)/koff2 + 
                  (c1*c2^2*kon^3)/(koff*koff2^2)))/kf2 - (1 + (c1*kon)/koff)*
                Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                   (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                     koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)] + 
               (koff2*(kf2*(1 + (c1*kon)/koff + (c2*kon)/koff2 - 
                    (c1*c2*kon^2)/(koff*koff2)) + (-1 + (c1*kon*(-1 + 
                       (c2*kon)/koff2))/koff)*Sqrt[koff^2*(1 + (c1*kon)/koff)^
                       2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                      (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                        koff)]))/kf2)))/(kf2*koff2*(kf2*(1 + koff2/kf2)*
              (-1 + (c2*kon)/kf2) + koff*(1 - (c2*kon)/kf2 + (c1*kon)/koff))*
            Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                (c1*kon*(-1 + (c2*kon)/koff2))/koff)]*(koff2 + c2*kon + 
             koff*(1 + (c1*kon)/koff) - Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
               koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*
            (-koff2 - c2*kon + kf*(2 + koff/kf - (c1*kon)/kf) + 
             Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                 2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                 (c1*kon*(-1 + (c2*kon)/koff2))/koff)])) + 
          ((1 + (c1*kon)/koff + (c2*kon)/koff2)*
            (koff^2*(1 + (c1*kon)/koff)^2 + kf2*(1 + koff2/kf2)*
              (koff2*(1 + (c2*kon)/koff2) - Sqrt[koff^2*(1 + (c1*kon)/koff)^
                   2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                  (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                    koff)]) + koff*(kf2*(-1 - (c1*kon)/koff + 
                 (koff2*(-2 - (c2*kon)/koff2 + (c1*kon*(-2 + (c2*kon)/koff2))/
                     koff))/kf2) + (1 + (c1*kon)/koff)*Sqrt[
                 koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                    2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + (c1*kon*
                      (-1 + (c2*kon)/koff2))/koff)]))*
            (c1^2*koff^2*kon^2*(1 + (c1*kon)/koff)^2 - c1*koff^2*kon*
              (kf*(1 + (c1*kon)/koff)^2 + koff2*(1 + (c1^3*kon^3)/koff^3 + 
                 (c1^2*kon^2*(4 - (3*c2*kon)/koff2))/koff^2 - 
                 (c1*kon*(-4 + (c2*kon)/koff2))/koff) + (c1*kon*
                 (1 + (c1*kon)/koff)*(kf2 + Sqrt[koff^2*(1 + (c1*kon)/koff)^
                      2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                     (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                       koff)]))/koff) + kf2*koff2*(-((c1*c2^2*koff2*kon^3*
                  (1 + (c2*kon)/koff2))/(kf2*koff)) + kf*(1 + (c1*kon)/koff)*
                Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                   (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                     koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)] + 
               (koff2^2*(kf*(1 + (c1*kon)/koff + (c2*kon)/koff2 + 
                    (c1*c2^2*kon^3)/(koff*koff2^2)) + (c1*c2*kon^2*
                    (kf2 + (c2*kf2*kon)/koff2 + (c2*kon*Sqrt[koff^2*
                          (1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/
                          koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                          (c1*kon*(-1 + (c2*kon)/koff2))/koff)])/koff2))/
                   (koff*koff2)))/kf2 + (koff2*(-((c1*c2*kf2*kon^2*
                     Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                        (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                          koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])/
                    (koff*koff2)) + kf*(kf2*(1 + (c1*kon)/koff + (c2*kon)/
                       koff2 - (c1*c2*kon^2)/(koff*koff2)) + 
                    (1 + (c1*kon)/koff - (c1*c2*kon^2)/(koff*koff2))*
                     Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                        (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                          koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])))/
                kf2) - kf2*koff*((c1*koff2^3*kon*(1 + (c1*kon)/koff - 
                  (c2*kon)/koff2 - (2*c2^2*kon^2)/koff2^2 - (c1*c2*kon^2)/
                   (koff*koff2) - (c2^3*kon^3)/koff2^3 + (3*c1*c2^2*kon^3)/
                   (koff*koff2^2)))/(kf2*koff) + (c1*kf*kon*Sqrt[
                  koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                     2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + (c1*kon*
                       (-1 + (c2*kon)/koff2))/koff)])/koff + (koff2^2*
                 (kf*(2 + (c2*kon)/koff2 - (2*c1^2*kon^2*(-1 + (c2*kon)/
                        koff2))/koff^2 + (c1*kon*(4 + (c2^2*kon^2)/koff2^2))/
                     koff) + (c1*kon*(kf2*(2 + (c1*kon)/koff + (2*c2*kon)/
                        koff2 + (c2^2*kon^2)/koff2^2 - (2*c1*c2*kon^2)/
                        (koff*koff2)) + (1 + (c1*kon)/koff + (c2^2*kon^2)/
                        koff2^2 - (2*c1*c2*kon^2)/(koff*koff2))*Sqrt[
                       koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/
                          koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                          (c1*kon*(-1 + (c2*kon)/koff2))/koff)]))/koff))/
                kf2 + (koff2*((c1*kf2*kon*((c1*kon)/koff - (c2*kon)/koff2)*
                    Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                       (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                         koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])/
                   koff + kf*(kf2*(1 + (c1*kon)/koff + (c1^2*kon^2)/koff^2 - 
                      (c1*c2*kon^2)/(koff*koff2)) + (1 + (c1^2*kon^2)/
                       koff^2 - (c1*kon*(-2 + (c2*kon)/koff2))/koff)*
                     Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                        (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                          koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])))/
                kf2) + koff^2*(-((c1*koff2^2*kon*(-2 + (c2*kon)/koff2 + 
                   (c2^2*kon^2)/koff2^2 + (c1^2*kon^2*(-2 + (3*c2*kon)/
                       koff2))/koff^2 + (c1*kon*(-4 + (2*c2*kon)/koff2 - 
                      (3*c2^2*kon^2)/koff2^2))/koff))/koff) + (c1*kon*
                 ((c1*kf2*kon*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                       (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                         koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])/
                   koff + kf*(1 + (c1*kon)/koff)*(kf2 + Sqrt[
                     koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/
                          koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                        (c1*kon*(-1 + (c2*kon)/koff2))/koff)])))/koff + koff2*
                (kf*(1 + (4*c1*kon)/koff + (c1^3*kon^3)/koff^3 - 
                   (2*c1^2*kon^2*(-2 + (c2*kon)/koff2))/koff^2) + 
                 (c1*kon*(kf2*(2 + (c2*kon)/koff2 + (c1^2*kon^2)/koff^2 - 
                      (2*c1*kon*(-1 + (c2*kon)/koff2))/koff) + 
                    (1 + (c1^2*kon^2)/koff^2 - (2*c1*kon*(-1 + (c2*kon)/
                          koff2))/koff)*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
                       koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                        (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                          koff)]))/koff))))/(kf2*koff2*
            Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                (c1*kon*(-1 + (c2*kon)/koff2))/koff)]*(koff2 + c2*kon + 
             koff*(1 + (c1*kon)/koff) + Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
               koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*
            ((c1*koff^2*kon*(1 - (c2*kon)/kf2 + (c1*kon)/koff))/kf + 
             kf2*(1 + koff2/kf2)*(-1 + (c2*kon)/kf2)*(-kf + koff2*
                (1 + (c2*kon)/koff2)) + (koff^2*(kf*(-1 + (c1*kon)/koff)*
                 (1 - (c2*kon)/kf2 + (c1*kon)/koff) + (c1*kf2*kon*
                  (-1 + (2*c2*koff2*kon)/kf2^2 + (koff2*(-2 - (c1*kon)/koff + 
                      (c2*kon)/koff2))/kf2))/koff))/kf + 
             (koff*(-((c1*kf2*koff2*(1 + koff2/kf2)*kon*(-1 + (c2*kon)/kf2))/
                  koff) - kf^2*(1 - (c2*kon)/kf2 + (c1*kon)/koff) + 
                kf*kf2*(1 - (c1*kon)/koff - (c2*koff2*kon*(2 - (c1*kon)/
                      koff + (c2*kon)/koff2))/kf2^2 + (2*koff2*(1 + 
                     (c1*c2*kon^2)/(koff*koff2)))/kf2)))/kf)))*kp*t)/
       (2*(1 + (c1*kon)/koff + (c2*kon)/koff2)^3), 
      (t*((2*c2*kon*((c1*c2*koff2*(1 + koff2/kf2)*kon^2)/koff - 
            c1*koff*kon*(koff2/kf2 - (c2*kon)/koff2 + (c1*koff2*kon)/(kf2*
                koff)) + koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2 + 
              (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/koff2^2 + 
              (koff2*(1 + (c1*kon)/koff + (c1^2*kon^2)/koff^2 - (c1*c2*kon^2)/
                  (koff*koff2)))/kf2)))/(koff*koff2*(1 + koff2/kf2)^2) - 
         (8*c1*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)*
           ((c1*c2*koff^2*kon^2)/(kf*koff2) + (c1*c2*koff2*kon^2)/koff + 
            (koff*koff2*(-(c2*kon*(1 + (c2*kon)/koff2)) + kf*
                (1 + (c1*kon)/koff + (c2*kon)/koff2 + (c1^2*kon^2)/koff^2 + 
                 (c2^2*kon^2)/koff2^2)))/kf + (koff^2*((c1*c2*kf*kon^2)/
                (koff*koff2) + koff2*(1 + (c2*kon)/koff2 + (c2^2*kon^2)/
                  koff2^2 - (c1*c2*kon^2)/(koff*koff2))))/kf)*
           (-koff^2 + 2*koff*koff2 - koff2^2 - 2*c1*koff*kon + 
            2*c2*koff*kon + 2*c1*koff2*kon - 2*c2*koff2*kon - c1^2*kon^2 - 
            2*c1*c2*kon^2 - c2^2*kon^2 - 2*kf*Sqrt[-4*koff*koff2*
                (1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] - 
            koff*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            koff2*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            c1*kon*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            c2*kon*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2]))/
          ((1 + koff/kf)^2*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/
                koff2) + (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/
                  koff2))^2]*(koff + koff2 + c1*kon + c2*kon - 
            Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
              (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2])*
           (2*kf + koff - koff2 - c1*kon - c2*kon + 
            Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
              (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2])*
           (koff + koff2 + c1*kon + c2*kon + Sqrt[-4*koff*koff2*(1 + 
                (c1*kon)/koff + (c2*kon)/koff2) + (koff*(1 + (c1*kon)/koff) + 
                koff2*(1 + (c2*kon)/koff2))^2])) + 
         (2*(c1*kon + c2*kon)*((c1*c2*koff2*kon^2)/koff + 
            (c2*koff^2*kon*((c1*kf*kon)/koff - (koff2*(kf2 - (c1*kf*kon)/
                   koff + (c1*kf2*kon)/koff))/kf2))/(kf*koff2) - 
            (koff*koff2*(-((c1*c2*kon^2)/koff) + kf*((c2*kon)/koff2 + 
                 (c1*kon*(1 + koff2/kf2 + (c2*kon)/kf2 + (2*c2*kon)/koff2))/
                  koff)))/kf))/(koff*(1 + koff/kf)*koff2*(1 + koff2/kf2)) + 
         (c2*kon*(c1*kon + c2*kon)*(1 + (c1*kon)/koff + (c2*kon)/koff2)*t)/
          (1 + koff2/kf2) - ((c1*kon + c2*kon)*(1 + (c1*kon)/koff + 
            (c2*kon)/koff2)*(c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*
                 kon)/koff))/kf)*t)/((1 + koff/kf)*(1 + koff2/kf2)) - 
         (4*c1*koff*koff2*kon*(c1*kon + c2*kon)*(1 + (c1*kon)/koff + 
             (c2*kon)/koff2)^2*(-koff^2 + 2*koff*koff2 - koff2^2 - 
            2*c1*koff*kon + 2*c2*koff*kon + 2*c1*koff2*kon - 2*c2*koff2*kon - 
            c1^2*kon^2 - 2*c1*c2*kon^2 - c2^2*kon^2 - 
            2*kf*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] - 
            koff*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            koff2*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            c1*kon*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            c2*kon*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2])*t)/
          ((1 + koff/kf)*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/
                koff2) + (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/
                  koff2))^2]*(koff + koff2 + c1*kon + c2*kon - 
            Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
              (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2])*
           (2*kf + koff - koff2 - c1*kon - c2*kon + 
            Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
              (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2])*
           (koff + koff2 + c1*kon + c2*kon + Sqrt[-4*koff*koff2*(1 + 
                (c1*kon)/koff + (c2*kon)/koff2) + (koff*(1 + (c1*kon)/koff) + 
                koff2*(1 + (c2*kon)/koff2))^2]))))/
       (2*(1 + (c1*kon)/koff + (c2*kon)/koff2)^3)}, 
     {-((((2*c1*c2*koff2*(1 + koff2/kf2)*kon^2*(1 + koff2/kf2 + 
             (c2*kon)/kf2))/(kf2*koff) + (c2*koff^4*kon*
            ((2*c1*kon*(1 + (c1*kon)/koff))/koff + 
             (koff2^2*(-1 + (c1*kon)/koff + (3*c2*kon)/koff2 + (2*c1^2*kon^2)/
                 koff^2 + (2*c2^2*kon^2)/koff2^2 + (4*c1*c2*kon^2)/
                 (koff*koff2)))/kf2^2 + (koff2*(-1 + (c2*kon)/koff2 + 
                (4*c1^2*kon^2)/koff^2 + (c1*kon*(3 + (2*c2*kon)/koff2))/
                 koff))/kf2))/(kf^3*koff2) + 
          (koff*koff2*((c2*kf*kon*(-1 - (c1*kon)/koff + (c2*kon)/koff2 + 
                (2*c1*koff2^2*kon*(1 + (c1*kon)/koff))/(kf2^2*koff) + 
                (koff2*(-1 + (c1*kon)/koff + (3*c2*kon)/koff2 + 
                   (2*c1^2*kon^2)/koff^2 + (2*c2^2*kon^2)/koff2^2 + 
                   (2*c1*c2*kon^2)/(koff*koff2)))/kf2))/koff2 + 
             (c1*kf2*(1 + koff2/kf2)*kon*(-1 + (c1*kon)/koff - (c2*kon)/
                 koff2 + (2*koff2*(-1 + (c1*kon)/koff + (c2*kon)/koff2 - 
                   (c2^2*kon^2)/koff2^2 + (c1*c2*kon^2)/(koff*koff2)))/kf2 + 
                (koff2^2*(-1 + (c1*kon)/koff + (3*c2*kon)/koff2 + 
                   (4*c2^2*kon^2)/koff2^2 + (2*c1*c2*kon^2)/(koff*koff2)))/
                 kf2^2))/koff))/(kf*kf2) + 
          (koff^2*((2*c1*c2*kf*kon^2)/(koff*koff2) + (c1*koff2^4*kon*(-1 + 
                (c2*kon)/koff2 + (2*c1^2*kon^2)/koff^2 + (2*c2^2*kon^2)/
                 koff2^2 + (c1*kon*(3 + (4*c2*kon)/koff2))/koff))/
              (kf2^3*koff) + (koff2*((c2*kf*kon*(-3 + (c1*kon)/koff + 
                   (3*c2*kon)/koff2 - (2*c1^2*kon^2)/koff^2 + (2*c1*c2*kon^2)/
                    (koff*koff2)))/koff2 + (c1*kf2*kon*(-1 + (3*c1*kon)/
                    koff + (c2*kon)/koff2 + (2*c1^2*kon^2)/koff^2 + 
                   (2*c2^2*kon^2)/koff2^2 + (2*c1*c2*kon^2)/(koff*koff2)))/
                 koff))/kf2 + (koff2^2*((c1*kf2*kon*(-3 + (9*c1*kon)/koff + 
                   (3*c2*kon)/koff2 + (6*c1^2*kon^2)/koff^2 + (8*c1*c2*kon^2)/
                    (koff*koff2)))/koff + (c2*kf*kon*(-3 + (3*c1*kon)/koff + 
                   (9*c2*kon)/koff2 + (6*c2^2*kon^2)/koff2^2 + 
                   (8*c1*c2*kon^2)/(koff*koff2)))/koff2))/kf2^2 + 
             (c1*koff2^3*kon*((2*c2*kf*kon*(2 + (c1*kon)/koff))/koff2 + 
                kf2*(-3 + (9*c1*kon)/koff + (3*c2*kon)/koff2 + (6*c1^2*kon^2)/
                   koff^2 + (10*c1*c2*kon^2)/(koff*koff2))))/(kf2^3*koff)))/
           kf^2 + (c2*koff^3*kon*((2*c1*(kf + kf2)*koff2^3*kon)/
              (kf2^3*koff) + (2*c1*kf*kon*(2 + (c1*kon)/koff))/koff + 
             (koff2*((2*c1*kf2*kon*(1 + (c2*kon)/koff2))/koff + 
                kf*(-3 + (5*c1*kon)/koff + (3*c2*kon)/koff2 + (2*c1^2*kon^2)/
                   koff^2 + (4*c1*c2*kon^2)/(koff*koff2))))/kf2 + 
             (koff2^2*((2*c1*kf2*kon*(2 + (c2*kon)/koff2))/koff + 
                kf*(-3 + (3*c1*kon)/koff + (9*c2*kon)/koff2 + (6*c2^2*kon^2)/
                   koff2^2 + (10*c1*c2*kon^2)/(koff*koff2))))/kf2^2))/
           (kf^3*koff2))*kp^2*t)/(koff*(1 + koff/kf)^3*koff2*
         (1 + koff2/kf2)^3*(1 + (c1*kon)/koff + (c2*kon)/koff2)^3)), 
      -(((-2*c2*kon)/(kf2*(1 + koff2/kf2)^2) - (2*c2*kon)/
           (koff2*(1 + koff2/kf2)^2) - (2*c1*c2*kon^2)/
           (kf2*koff*(1 + koff2/kf2)^2) + (2*c1*c2*kon^2)/
           (koff2^2*(1 + koff2/kf2)^2) + (4*c1*c2*kon^2)/
           (kf2*koff2*(1 + koff2/kf2)^2) + (2*c2^2*kon^2)/
           (kf2*koff2*(1 + koff2/kf2)^2) - (2*c1*c2*kon^2)/
           (koff*koff2*(1 + koff2/kf2)^2) - (2*c1^2*c2*kon^3)/
           (kf2*koff^2*(1 + koff2/kf2)^2) + (2*c1*c2^2*kon^3)/
           (kf2*koff2^2*(1 + koff2/kf2)^2) + (2*c2^3*kon^3)/
           (kf2*koff2^2*(1 + koff2/kf2)^2) + (2*c1^2*c2*kon^3)/
           (koff*koff2^2*(1 + koff2/kf2)^2) + (2*c1*c2^2*kon^3)/
           (koff*koff2^2*(1 + koff2/kf2)^2) - (2*c1^2*c2*kon^3)/
           (koff^2*koff2*(1 + koff2/kf2)^2) + (4*c1^2*c2*kon^3)/
           (kf2*koff*koff2*(1 + koff2/kf2)^2) + (4*c1*c2^2*kon^3)/
           (kf2*koff*koff2*(1 + koff2/kf2)^2) - (2*c1*c2^2*kon^3)/
           (koff*koff2*(koff + (koff*koff2)/kf2)) + 
          (2*c1*(kf + koff - koff2)*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)^
             2*(c2*kon + kf*(-1 + (koff*(-1 + (c1*kon)/koff))/kf)))/
           (kf^2*(1 + koff/kf)^2*(kf*(1 + koff/kf)*(-1 + (c1*kon)/kf) + 
             koff2*(1 - (c1*kon)/kf + (c2*kon)/koff2))) + 
          (2*(c1*kon + c2*kon)*((c2*koff*kon*(1 + koff/kf + (c1*kon)/kf))/
              koff2 + (c1*koff2^2*kon*(1 + (c2*kon)/koff2))/(kf2*koff) - 
             (c1*koff2*kon*((c2*koff*kon)/koff2 + kf2*(-1 + (c2*koff*kon)/
                   (kf*koff2))))/(kf2*koff)))/(koff*(1 + koff/kf)*koff2*
            (1 + koff2/kf2)) + (2*kf*(1 + (c1*kon)/koff + (c2*kon)/koff2)*
            (koff^2*(1 + (c1*kon)/koff)^2 + kf2*(1 + koff2/kf2)*
              (koff2*(1 + (c2*kon)/koff2) + Sqrt[koff^2*(1 + (c1*kon)/koff)^
                   2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                  (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                    koff)]) - koff*(kf2*(1 + (c1*kon)/koff + 
                 (koff2*(2 + (c2*kon)/koff2 - (c1*kon*(-2 + (c2*kon)/koff2))/
                     koff))/kf2) + (1 + (c1*kon)/koff)*Sqrt[
                 koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                    2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + (c1*kon*
                      (-1 + (c2*kon)/koff2))/koff)]))*
            (-(c1*koff^2*kon*(1 + (c1*kon)/koff)^2) + 
             koff^2*(koff2*(1 + (4*c1*kon)/koff + (c1^3*kon^3)/koff^3 - 
                 (2*c1^2*kon^2*(-2 + (c2*kon)/koff2))/koff^2) - (c1*kon*
                 (1 + (c1*kon)/koff)*(-kf2 + Sqrt[koff^2*(1 + (c1*kon)/koff)^
                      2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                     (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                       koff)]))/koff) - kf2*koff*((koff2^2*(2 + (c2*kon)/
                   koff2 - (2*c1^2*kon^2*(-1 + (c2*kon)/koff2))/koff^2 + 
                  (c1*kon*(4 + (c2^2*kon^2)/koff2^2))/koff))/kf2 - 
               (c1*kon*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                    (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                      koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])/koff + 
               (koff2*(kf2*(1 + (c1*kon)/koff + (c1^2*kon^2)/koff^2 - 
                    (c1*c2*kon^2)/(koff*koff2)) - (1 + (c1^2*kon^2)/koff^2 - 
                    (c1*kon*(-2 + (c2*kon)/koff2))/koff)*Sqrt[
                    koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/
                         koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                       (c1*kon*(-1 + (c2*kon)/koff2))/koff)]))/kf2) + 
             kf2*koff2*((koff2^2*(1 + (c1*kon)/koff + (c2*kon)/koff2 + 
                  (c1*c2^2*kon^3)/(koff*koff2^2)))/kf2 - (1 + (c1*kon)/koff)*
                Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                   (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                     koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)] + 
               (koff2*(kf2*(1 + (c1*kon)/koff + (c2*kon)/koff2 - 
                    (c1*c2*kon^2)/(koff*koff2)) + (-1 + (c1*kon*(-1 + 
                       (c2*kon)/koff2))/koff)*Sqrt[koff^2*(1 + (c1*kon)/koff)^
                       2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                      (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                        koff)]))/kf2)))/(kf2*koff2*(kf2*(1 + koff2/kf2)*
              (-1 + (c2*kon)/kf2) + koff*(1 - (c2*kon)/kf2 + (c1*kon)/koff))*
            Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                (c1*kon*(-1 + (c2*kon)/koff2))/koff)]*(koff2 + c2*kon + 
             koff*(1 + (c1*kon)/koff) - Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
               koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*
            (-koff2 - c2*kon + kf*(2 + koff/kf - (c1*kon)/kf) + 
             Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                 2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                 (c1*kon*(-1 + (c2*kon)/koff2))/koff)])) + 
          ((1 + (c1*kon)/koff + (c2*kon)/koff2)*
            (koff^2*(1 + (c1*kon)/koff)^2 + kf2*(1 + koff2/kf2)*
              (koff2*(1 + (c2*kon)/koff2) - Sqrt[koff^2*(1 + (c1*kon)/koff)^
                   2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                  (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                    koff)]) + koff*(kf2*(-1 - (c1*kon)/koff + 
                 (koff2*(-2 - (c2*kon)/koff2 + (c1*kon*(-2 + (c2*kon)/koff2))/
                     koff))/kf2) + (1 + (c1*kon)/koff)*Sqrt[
                 koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                    2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + (c1*kon*
                      (-1 + (c2*kon)/koff2))/koff)]))*
            (c1^2*koff^2*kon^2*(1 + (c1*kon)/koff)^2 - c1*koff^2*kon*
              (kf*(1 + (c1*kon)/koff)^2 + koff2*(1 + (c1^3*kon^3)/koff^3 + 
                 (c1^2*kon^2*(4 - (3*c2*kon)/koff2))/koff^2 - 
                 (c1*kon*(-4 + (c2*kon)/koff2))/koff) + (c1*kon*
                 (1 + (c1*kon)/koff)*(kf2 + Sqrt[koff^2*(1 + (c1*kon)/koff)^
                      2 + koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                     (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                       koff)]))/koff) + kf2*koff2*(-((c1*c2^2*koff2*kon^3*
                  (1 + (c2*kon)/koff2))/(kf2*koff)) + kf*(1 + (c1*kon)/koff)*
                Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                   (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                     koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)] + 
               (koff2^2*(kf*(1 + (c1*kon)/koff + (c2*kon)/koff2 + 
                    (c1*c2^2*kon^3)/(koff*koff2^2)) + (c1*c2*kon^2*
                    (kf2 + (c2*kf2*kon)/koff2 + (c2*kon*Sqrt[koff^2*
                          (1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/
                          koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                          (c1*kon*(-1 + (c2*kon)/koff2))/koff)])/koff2))/
                   (koff*koff2)))/kf2 + (koff2*(-((c1*c2*kf2*kon^2*
                     Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                        (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                          koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])/
                    (koff*koff2)) + kf*(kf2*(1 + (c1*kon)/koff + (c2*kon)/
                       koff2 - (c1*c2*kon^2)/(koff*koff2)) + 
                    (1 + (c1*kon)/koff - (c1*c2*kon^2)/(koff*koff2))*
                     Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                        (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                          koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])))/
                kf2) - kf2*koff*((c1*koff2^3*kon*(1 + (c1*kon)/koff - 
                  (c2*kon)/koff2 - (2*c2^2*kon^2)/koff2^2 - (c1*c2*kon^2)/
                   (koff*koff2) - (c2^3*kon^3)/koff2^3 + (3*c1*c2^2*kon^3)/
                   (koff*koff2^2)))/(kf2*koff) + (c1*kf*kon*Sqrt[
                  koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                     2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + (c1*kon*
                       (-1 + (c2*kon)/koff2))/koff)])/koff + (koff2^2*
                 (kf*(2 + (c2*kon)/koff2 - (2*c1^2*kon^2*(-1 + (c2*kon)/
                        koff2))/koff^2 + (c1*kon*(4 + (c2^2*kon^2)/koff2^2))/
                     koff) + (c1*kon*(kf2*(2 + (c1*kon)/koff + (2*c2*kon)/
                        koff2 + (c2^2*kon^2)/koff2^2 - (2*c1*c2*kon^2)/
                        (koff*koff2)) + (1 + (c1*kon)/koff + (c2^2*kon^2)/
                        koff2^2 - (2*c1*c2*kon^2)/(koff*koff2))*Sqrt[
                       koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/
                          koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                          (c1*kon*(-1 + (c2*kon)/koff2))/koff)]))/koff))/
                kf2 + (koff2*((c1*kf2*kon*((c1*kon)/koff - (c2*kon)/koff2)*
                    Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                       (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                         koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])/
                   koff + kf*(kf2*(1 + (c1*kon)/koff + (c1^2*kon^2)/koff^2 - 
                      (c1*c2*kon^2)/(koff*koff2)) + (1 + (c1^2*kon^2)/
                       koff^2 - (c1*kon*(-2 + (c2*kon)/koff2))/koff)*
                     Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                        (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                          koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])))/
                kf2) + koff^2*(-((c1*koff2^2*kon*(-2 + (c2*kon)/koff2 + 
                   (c2^2*kon^2)/koff2^2 + (c1^2*kon^2*(-2 + (3*c2*kon)/
                       koff2))/koff^2 + (c1*kon*(-4 + (2*c2*kon)/koff2 - 
                      (3*c2^2*kon^2)/koff2^2))/koff))/koff) + (c1*kon*
                 ((c1*kf2*kon*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*
                       (1 + (c2*kon)/koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/
                         koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])/
                   koff + kf*(1 + (c1*kon)/koff)*(kf2 + Sqrt[
                     koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/
                          koff2)^2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                        (c1*kon*(-1 + (c2*kon)/koff2))/koff)])))/koff + koff2*
                (kf*(1 + (4*c1*kon)/koff + (c1^3*kon^3)/koff^3 - 
                   (2*c1^2*kon^2*(-2 + (c2*kon)/koff2))/koff^2) + 
                 (c1*kon*(kf2*(2 + (c2*kon)/koff2 + (c1^2*kon^2)/koff^2 - 
                      (2*c1*kon*(-1 + (c2*kon)/koff2))/koff) + 
                    (1 + (c1^2*kon^2)/koff^2 - (2*c1*kon*(-1 + (c2*kon)/
                          koff2))/koff)*Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
                       koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                        (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/
                          koff)]))/koff))))/(kf2*koff2*
            Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + koff2^2*(1 + (c2*kon)/koff2)^
                2 + 2*koff*koff2*(-1 - (c2*kon)/koff2 + 
                (c1*kon*(-1 + (c2*kon)/koff2))/koff)]*(koff2 + c2*kon + 
             koff*(1 + (c1*kon)/koff) + Sqrt[koff^2*(1 + (c1*kon)/koff)^2 + 
               koff2^2*(1 + (c2*kon)/koff2)^2 + 2*koff*koff2*
                (-1 - (c2*kon)/koff2 + (c1*kon*(-1 + (c2*kon)/koff2))/koff)])*
            ((c1*koff^2*kon*(1 - (c2*kon)/kf2 + (c1*kon)/koff))/kf + 
             kf2*(1 + koff2/kf2)*(-1 + (c2*kon)/kf2)*(-kf + koff2*
                (1 + (c2*kon)/koff2)) + (koff^2*(kf*(-1 + (c1*kon)/koff)*
                 (1 - (c2*kon)/kf2 + (c1*kon)/koff) + (c1*kf2*kon*
                  (-1 + (2*c2*koff2*kon)/kf2^2 + (koff2*(-2 - (c1*kon)/koff + 
                      (c2*kon)/koff2))/kf2))/koff))/kf + 
             (koff*(-((c1*kf2*koff2*(1 + koff2/kf2)*kon*(-1 + (c2*kon)/kf2))/
                  koff) - kf^2*(1 - (c2*kon)/kf2 + (c1*kon)/koff) + 
                kf*kf2*(1 - (c1*kon)/koff - (c2*koff2*kon*(2 - (c1*kon)/
                      koff + (c2*kon)/koff2))/kf2^2 + (2*koff2*(1 + 
                     (c1*c2*kon^2)/(koff*koff2)))/kf2)))/kf)))*kp*t)/
       (2*(1 + (c1*kon)/koff + (c2*kon)/koff2)^3), 
      (kp*((c1*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)^2)/
          (koff*(1 + koff/kf)) + (c2*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)^
            2)/(koff2*(1 + koff2/kf2)) + 
         (2*((c2*koff*(1 + koff/kf)*kon*(1 + koff/kf + (c1*kon)/kf)^2)/
             koff2 + (c1*koff2^4*kon*((1 + (c2*kon)/koff2)^2 + (koff^2*
                 (1 + (c1*kon)/koff + (c2*kon)/koff2)^2)/kf^2 + 
               (koff*(1 + (c2*kon)/koff2)*(2 + (c1*kon)/koff + (2*c2*kon)/
                   koff2))/kf))/(kf2^3*koff) + (c1*koff2^3*kon*
              ((c1*c2*(1 + koff/kf)*kon^2)/koff2 + kf2*(3 + (4*c2*kon)/
                  koff2 + (c2^2*kon^2)/koff2^2 + (koff*(6 + (3*c1*kon)/koff + 
                    (7*c2*kon)/koff2 + (c2^2*kon^2)/koff2^2 + (2*c1*c2*kon^2)/
                     (koff*koff2)))/kf + (koff^2*(3 + (6*c1*kon)/koff + 
                    (3*c2*kon)/koff2 + (3*c1^2*kon^2)/koff^2 + 
                    (5*c1*c2*kon^2)/(koff*koff2)))/kf^2)))/(kf2^3*koff) + 
            (koff2*((c2*koff*(1 + koff/kf)*kon*(1 + koff/kf + (c1*kon)/kf)*
                 (2 - (c1*kon)/koff + (c2*kon)/koff2 + (koff*(2 + (2*c1*kon)/
                      koff + (c2*kon)/koff2))/kf))/koff2 + (c1*kf2*kon*
                 (1 + (c2^2*koff^3*kon^2)/(kf^3*koff2^2) + 
                  (koff*(2 + (c1*kon)/koff - (c2*kon)/koff2))/kf + 
                  (koff^2*(1 - (c2*kon)/koff2 + (c1^2*kon^2)/koff^2 + 
                     (c2^2*kon^2)/koff2^2 + (c1*kon*(2 + (c2*kon)/koff2))/
                      koff))/kf^2))/koff))/kf2 + 
            (koff2^2*((c2*koff*(1 + koff/kf)*kon*((c1^2*kon^2)/koff^2 + 
                  (c1*kon*(-1 + (c2*kon)/koff2))/koff + (1 + (c2*kon)/koff2)^
                   2 + (koff^2*(1 + (c1*kon)/koff + (c2*kon)/koff2)^2)/kf^2 + 
                  (koff*((c1*kon)/koff - (c1^2*kon^2)/koff^2 + 
                     (3*c1*c2*kon^2)/(koff*koff2) + 2*(1 + (c2*kon)/koff2)^
                       2))/kf))/koff2 + (c1*kf2*kon*(3 + (2*c2*kon)/koff2 + 
                  (c2^2*koff^3*kon^2)/(kf^3*koff2^2) + (koff*(6 + (2*c2*kon)/
                      koff2 - (c2^2*kon^2)/koff2^2 + (c1*kon*(3 + (c2*kon)/
                         koff2))/koff))/kf + (koff^2*(3 + (3*c1^2*kon^2)/
                      koff^2 + (c1*kon*(6 + (4*c2*kon)/koff2))/koff))/kf^2))/
                koff))/kf2^2)*kp)/(koff*(1 + koff/kf)^3*koff2*
           (1 + koff2/kf2)^3))*t)/(1 + (c1*kon)/koff + (c2*kon)/koff2)^3, 
      (kp*t*((-2*c2*(1 + koff/kf)^2*kon*((c1*koff2*(1 + koff2/kf2)*kon*
              (1 + koff2/kf2 + (c2*kon)/kf2))/koff + 
            (c2*koff^2*kon*(1 + (c1*kon)/koff + (koff2*(2 + (2*c1*kon)/koff + 
                  (c2*kon)/koff2))/kf2))/(kf*koff2) + 
            (koff*((c2*kf*kon)/koff2 + (c1*koff2^2*kon*(kf + (c1*kf*kon)/
                   koff - (c2*kf2*kon)/koff2))/(kf2^2*koff) + (koff2*
                 (-((c1*c2*kf2*kon^2)/(koff*koff2)) + kf*((c1*kon)/koff + 
                    (2*c2*kon)/koff2 + (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/
                     koff2^2 + (c1*c2*kon^2)/(koff*koff2))))/kf2))/kf))/
          (koff*koff2) + (2*c2*(1 + koff/kf)^2*kon*
           ((c1*c2*koff2*(1 + koff2/kf2)*kon^2)/koff + 
            (koff^2*(-((c1*kf*kon*(1 + (c1*kon)/koff))/koff) + (koff2^3*
                 (1 + (c1*kon)/koff + (c2*kon)/koff2)^2)/kf2^2 + (koff2^2*
                 (kf2 - (c1*kf*kon)/koff + (c1*kf2*kon)/koff)*
                 (2 + (2*c1*kon)/koff + (c2*kon)/koff2))/kf2^2 + 
               (koff2*(kf2*(1 + (c1*kon)/koff)^2 - (c1*kf*kon*(3 + (3*c1*kon)/
                      koff + (c2*kon)/koff2))/koff))/kf2))/kf + 
            (koff*koff2*((c1*c2*(1 + koff2/kf2)*kon^2)/koff + kf*
                (1 + (c1*kon)/koff + (c1^2*kon^2)/koff^2 - (c1*c2*kon^2)/
                  (koff*koff2) + (koff2*(2 + (c2*kon)/koff2 + (2*c1^2*kon^2)/
                     koff^2 - (c1*kon*(-2 + (c2*kon)/koff2))/koff))/kf2 + 
                 (koff2^2*((c1^2*kon^2)/koff^2 + (c1*kon*(1 + (c2*kon)/
                        koff2))/koff + (1 + (c2*kon)/koff2)^2))/kf2^2)))/kf))/
          (koff*koff2^2) - (2*c1*(1 + koff2/kf2)^2*kon*
           ((c2*koff*(1 + koff/kf)*kon*(1 + koff/kf + (c1*kon)/kf))/koff2 + 
            (c1*koff2^2*kon*(1 + (c2*kon)/koff2 + (koff*(2 + (c1*kon)/koff + 
                  (2*c2*kon)/koff2))/kf))/(kf2*koff) + 
            (koff2*(-((c1*c2*(1 + koff/kf)*kon^2)/koff2) + kf2*
                ((c1^2*kon^2)/(kf*koff) + (c2*koff*(1 + koff/kf)*kon*
                   (1 + (c2*kon)/koff2))/(kf*koff2) + (c1*kon*(1 + 
                    (koff*(2 + (c2*kon)/koff2))/kf))/koff)))/kf2))/
          (koff*koff2) + (2*c1*(1 + koff2/kf2)^2*kon*
           (-(c2*koff2*kon*(1 + (c2*kon)/koff2)) + 
            (koff*koff2*(-(c2*kon*(3 + (c1*kon)/koff + (3*c2*kon)/koff2)) + 
               kf*(1 + (c2*kon)/koff2 + (c2^2*kon^2)/koff2^2 - (c1*c2*kon^2)/
                  (koff*koff2) + (koff2*(1 + (c2*kon)/koff2)^2)/kf2)))/kf + 
            (koff^3*((c1*c2*kf*kon^2)/(koff*koff2) + (koff2^2*
                 (1 + (c1*kon)/koff + (c2*kon)/koff2)^2)/kf2 + (koff2*
                 ((c1*c2*kf*kon^2)/(koff*koff2) + kf2*(1 + (c2*kon)/koff2 + 
                    (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/koff2^2 + 
                    (c1*kon*(2 + (c2*kon)/koff2))/koff)))/kf2))/kf^2 + 
            (koff^2*((c1*c2*kf^2*kon^2)/(koff*koff2) + (koff2^2*
                 (2 + (c1*kon)/koff + (2*c2*kon)/koff2)*(kf + (c2*kf*kon)/
                   koff2 - (c2*kf2*kon)/koff2))/kf2 + (kf*koff2*
                 ((c1*c2*kf*kon^2)/(koff*koff2) + kf2*((c1*kon)/koff - 
                    (c1*c2*kon^2)/(koff*koff2) + 2*(1 + (c2*kon)/koff2 + 
                      (c2^2*kon^2)/koff2^2))))/kf2))/kf^2))/(koff^2*koff2) + 
         c2*(1 + koff/kf)^2*(1 + koff2/kf2)*kon*(1 + (c1*kon)/koff + 
           (c2*kon)/koff2)*((c1*kon)/koff + (c2*kon)/koff2 + 
           (c2*koff*kon)/(kf*koff2) + (c1*koff2*kon)/(kf2*koff))*t + 
         c1*(1 + koff/kf)*(1 + koff2/kf2)^2*kon*(1 + (c1*kon)/koff + 
           (c2*kon)/koff2)*((c1*kon)/koff + (c2*kon)/koff2 + 
           (c2*koff*kon)/(kf*koff2) + (c1*koff2*kon)/(kf2*koff))*t + 
         (c2*(1 + koff/kf)^2*(1 + koff2/kf2)*kon*(1 + (c1*kon)/koff + 
            (c2*kon)/koff2)*(c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*
                 kon)/koff))/kf)*t)/koff2 + (c1*(1 + koff/kf)*
           (1 + koff2/kf2)^2*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)*
           (c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*kon)/koff))/kf)*
           t)/koff - 2*(1 + koff/kf)*(1 + koff2/kf2)*(1 + (c1*kon)/koff + 
           (c2*kon)/koff2)*((c1*kon)/koff + (c2*kon)/koff2 + 
           (c2*koff*kon)/(kf*koff2) + (c1*koff2*kon)/(kf2*koff))*
          (c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*kon)/koff))/kf)*
          t))/(2*(1 + koff/kf)^3*(1 + koff2/kf2)^3*
        (1 + (c1*kon)/koff + (c2*kon)/koff2)^3)}, 
     {(kp*t*((2*c2*(1 + koff/kf)^2*kon*((c1*c2*koff2*(1 + koff2/kf2)*kon^2)/
             (kf2*koff) + (koff^2*((c1*c2*kon^2)/(koff*koff2) + (koff2^2*
                 (1 + (c1*kon)/koff))/kf2^2 + (koff2*(1 + (c1*kon)/koff + 
                  (c2*kon)/koff2 + (c2^2*kon^2)/koff2^2 + (2*c1*c2*kon^2)/
                   (koff*koff2)))/kf2))/kf + 
            (koff*koff2*(-((c1*kf2*(1 + koff2/kf2)*kon*(1 + koff2/kf2 + 
                   (c2*kon)/koff2))/koff) + kf*(1 + (c2*kon)/koff2 + 
                 (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/koff2^2 + 
                 (koff2*(1 + (c1*kon)/koff)^2)/kf2 + (c1*kon*(2 + (c2*kon)/
                     koff2))/koff)))/(kf*kf2)))/(koff*koff2) + 
         (2*c1*(1 + koff2/kf2)^2*kon*((c1*c2*koff2*kon^2)/(kf2*koff) + 
            (c2*koff^3*kon*(-(koff2/kf2) + (c1*kon)/koff))/(kf^2*koff2) + 
            (koff^2*((c1*c2*kf*kon^2)/(koff*koff2) + (koff2^2*
                 (1 + (c2*kon)/koff2))/kf2 + (koff2*(-((c2*kf*kon*
                     (2 + (c1*kon)/koff))/koff2) + kf2*(1 + (c2*kon)/koff2)^
                    2))/kf2))/kf^2 + (koff*koff2*(-((c2*kf*kon*
                  (1 + (c1*kon)/koff))/koff2) + kf2*((c1^2*kon^2)/koff^2 + 
                 (c1*kon*(1 + (c2*kon)/koff2))/koff + (1 + (c2*kon)/koff2)^
                  2 + (koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2 + 
                    (c1^2*kon^2)/koff^2 + (2*c1*c2*kon^2)/(koff*koff2)))/
                  kf2)))/(kf*kf2)))/(koff*koff2) - 
         (2*c1*(1 + koff2/kf2)^2*kon*(-((c1*c2*koff2*kon^2)/koff) + 
            (c2*koff^3*kon*(-((c1*kf*kon)/koff) + (koff2*(kf2 - (c1*kf*kon)/
                   koff + (c1*kf2*kon)/koff))/kf2))/(kf^2*koff2) + 
            (koff*koff2*(c2*kon*(1 - (c1*kon)/koff + (c2*kon)/koff2) + kf*
                ((c2*kon)/koff2 + (c1*kon*(1 + koff2/kf2 + (c2*kon)/kf2 + 
                    (2*c2*kon)/koff2))/koff)))/kf + 
            (koff^2*(-((c1*c2*kf^2*kon^2)/(koff*koff2)) + (koff2^2*
                 ((c2*kf2*kon*(1 + (c2*kon)/koff2))/koff2 + (c1*kf*kon*
                    (2 + (c1*kon)/koff + (2*c2*kon)/koff2))/koff))/kf2 + 
               (kf*koff2*(-((c1*c2*kf*kon^2)/(koff*koff2)) + 
                  kf2*((2*c1*kon)/koff + (2*c2*kon)/koff2 + (c1^2*kon^2)/
                     koff^2 + (4*c1*c2*kon^2)/(koff*koff2))))/kf2))/kf^2))/
          (kf*koff*koff2) - (2*c2*(1 + koff/kf)^2*kon*
           (-((c1*c2*koff2*(1 + koff2/kf2)*kon^2)/koff) + 
            (koff^2*(-((c1*c2*kf*kon^2)/(koff*koff2)) + (koff2*
                 ((c2*kf2*kon*(1 + (c1*kon)/koff))/koff2 + (c1*kf*kon*
                    (1 + (c1*kon)/koff - (c2*kon)/koff2))/koff))/kf2 + 
               (koff2^2*((c1*kf*kon*(1 + (c1*kon)/koff))/koff + 
                  (c2*kf2*kon*(2 + (2*c1*kon)/koff + (c2*kon)/koff2))/koff2))/
                kf2^2))/kf + (koff*koff2*(-((c1*c2*(1 + koff2/kf2)*kon^2)/
                 koff) + kf*((c2*kon*(1 + (koff2*(2 + (c2*kon)/koff2))/kf2))/
                  koff2 + (c1*kon*(1 + (2*c2*kon)/koff2 + (koff2^2*
                      (1 + (c2*kon)/koff2))/kf2^2 + (koff2*(2 + (4*c2*kon)/
                        koff2))/kf2))/koff)))/kf))/(kf2*koff*koff2) + 
         c2*(1 + koff/kf)^2*(1 + koff2/kf2)*kon*(1 + (c1*kon)/koff + 
           (c2*kon)/koff2)*((c2*kon)/kf2 + 
           (koff*((c2*kon)/kf2 + (c1*kon)/koff + (c1*koff2*kon)/(kf2*koff)))/
            kf)*t + c1*(1 + koff/kf)*(1 + koff2/kf2)^2*kon*
          (1 + (c1*kon)/koff + (c2*kon)/koff2)*((c2*kon)/kf2 + 
           (koff*((c2*kon)/kf2 + (c1*kon)/koff + (c1*koff2*kon)/(kf2*koff)))/
            kf)*t + (c2*(1 + koff/kf)^2*(1 + koff2/kf2)*kon*
           (1 + (c1*kon)/koff + (c2*kon)/koff2)*(c2*kon + 
            (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*kon)/koff))/kf)*t)/kf2 + 
         (c1*(1 + koff/kf)*(1 + koff2/kf2)^2*kon*(1 + (c1*kon)/koff + 
            (c2*kon)/koff2)*(c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*
                 kon)/koff))/kf)*t)/kf - 2*(1 + koff/kf)*(1 + koff2/kf2)*
          (1 + (c1*kon)/koff + (c2*kon)/koff2)*((c2*kon)/kf2 + 
           (koff*((c2*kon)/kf2 + (c1*kon)/koff + (c1*koff2*kon)/(kf2*koff)))/
            kf)*(c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*kon)/koff))/
            kf)*t))/(2*(1 + koff/kf)^3*(1 + koff2/kf2)^3*
        (1 + (c1*kon)/koff + (c2*kon)/koff2)^3), 
      (t*((2*c2*kon*((c1*c2*koff2*(1 + koff2/kf2)*kon^2)/koff - 
            c1*koff*kon*(koff2/kf2 - (c2*kon)/koff2 + (c1*koff2*kon)/(kf2*
                koff)) + koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2 + 
              (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/koff2^2 + 
              (koff2*(1 + (c1*kon)/koff + (c1^2*kon^2)/koff^2 - (c1*c2*kon^2)/
                  (koff*koff2)))/kf2)))/(koff*koff2*(1 + koff2/kf2)^2) - 
         (8*c1*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)*
           ((c1*c2*koff^2*kon^2)/(kf*koff2) + (c1*c2*koff2*kon^2)/koff + 
            (koff*koff2*(-(c2*kon*(1 + (c2*kon)/koff2)) + kf*
                (1 + (c1*kon)/koff + (c2*kon)/koff2 + (c1^2*kon^2)/koff^2 + 
                 (c2^2*kon^2)/koff2^2)))/kf + (koff^2*((c1*c2*kf*kon^2)/
                (koff*koff2) + koff2*(1 + (c2*kon)/koff2 + (c2^2*kon^2)/
                  koff2^2 - (c1*c2*kon^2)/(koff*koff2))))/kf)*
           (-koff^2 + 2*koff*koff2 - koff2^2 - 2*c1*koff*kon + 
            2*c2*koff*kon + 2*c1*koff2*kon - 2*c2*koff2*kon - c1^2*kon^2 - 
            2*c1*c2*kon^2 - c2^2*kon^2 - 2*kf*Sqrt[-4*koff*koff2*
                (1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] - 
            koff*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            koff2*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            c1*kon*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            c2*kon*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2]))/
          ((1 + koff/kf)^2*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/
                koff2) + (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/
                  koff2))^2]*(koff + koff2 + c1*kon + c2*kon - 
            Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
              (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2])*
           (2*kf + koff - koff2 - c1*kon - c2*kon + 
            Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
              (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2])*
           (koff + koff2 + c1*kon + c2*kon + Sqrt[-4*koff*koff2*(1 + 
                (c1*kon)/koff + (c2*kon)/koff2) + (koff*(1 + (c1*kon)/koff) + 
                koff2*(1 + (c2*kon)/koff2))^2])) + 
         (2*(c1*kon + c2*kon)*((c1*c2*koff2*kon^2)/koff + 
            (c2*koff^2*kon*((c1*kf*kon)/koff - (koff2*(kf2 - (c1*kf*kon)/
                   koff + (c1*kf2*kon)/koff))/kf2))/(kf*koff2) - 
            (koff*koff2*(-((c1*c2*kon^2)/koff) + kf*((c2*kon)/koff2 + 
                 (c1*kon*(1 + koff2/kf2 + (c2*kon)/kf2 + (2*c2*kon)/koff2))/
                  koff)))/kf))/(koff*(1 + koff/kf)*koff2*(1 + koff2/kf2)) + 
         (c2*kon*(c1*kon + c2*kon)*(1 + (c1*kon)/koff + (c2*kon)/koff2)*t)/
          (1 + koff2/kf2) - ((c1*kon + c2*kon)*(1 + (c1*kon)/koff + 
            (c2*kon)/koff2)*(c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*
                 kon)/koff))/kf)*t)/((1 + koff/kf)*(1 + koff2/kf2)) - 
         (4*c1*koff*koff2*kon*(c1*kon + c2*kon)*(1 + (c1*kon)/koff + 
             (c2*kon)/koff2)^2*(-koff^2 + 2*koff*koff2 - koff2^2 - 
            2*c1*koff*kon + 2*c2*koff*kon + 2*c1*koff2*kon - 2*c2*koff2*kon - 
            c1^2*kon^2 - 2*c1*c2*kon^2 - c2^2*kon^2 - 
            2*kf*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] - 
            koff*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            koff2*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            c1*kon*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2] + 
            c2*kon*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
               (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2])*t)/
          ((1 + koff/kf)*Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/
                koff2) + (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/
                  koff2))^2]*(koff + koff2 + c1*kon + c2*kon - 
            Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
              (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2])*
           (2*kf + koff - koff2 - c1*kon - c2*kon + 
            Sqrt[-4*koff*koff2*(1 + (c1*kon)/koff + (c2*kon)/koff2) + 
              (koff*(1 + (c1*kon)/koff) + koff2*(1 + (c2*kon)/koff2))^2])*
           (koff + koff2 + c1*kon + c2*kon + Sqrt[-4*koff*koff2*(1 + 
                (c1*kon)/koff + (c2*kon)/koff2) + (koff*(1 + (c1*kon)/koff) + 
                koff2*(1 + (c2*kon)/koff2))^2]))))/
       (2*(1 + (c1*kon)/koff + (c2*kon)/koff2)^3), 
      (kp*t*((-2*c2*(1 + koff/kf)^2*kon*((c1*koff2*(1 + koff2/kf2)*kon*
              (1 + koff2/kf2 + (c2*kon)/kf2))/koff + 
            (c2*koff^2*kon*(1 + (c1*kon)/koff + (koff2*(2 + (2*c1*kon)/koff + 
                  (c2*kon)/koff2))/kf2))/(kf*koff2) + 
            (koff*((c2*kf*kon)/koff2 + (c1*koff2^2*kon*(kf + (c1*kf*kon)/
                   koff - (c2*kf2*kon)/koff2))/(kf2^2*koff) + (koff2*
                 (-((c1*c2*kf2*kon^2)/(koff*koff2)) + kf*((c1*kon)/koff + 
                    (2*c2*kon)/koff2 + (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/
                     koff2^2 + (c1*c2*kon^2)/(koff*koff2))))/kf2))/kf))/
          (koff*koff2) + (2*c2*(1 + koff/kf)^2*kon*
           ((c1*c2*koff2*(1 + koff2/kf2)*kon^2)/koff + 
            (koff^2*(-((c1*kf*kon*(1 + (c1*kon)/koff))/koff) + (koff2^3*
                 (1 + (c1*kon)/koff + (c2*kon)/koff2)^2)/kf2^2 + (koff2^2*
                 (kf2 - (c1*kf*kon)/koff + (c1*kf2*kon)/koff)*
                 (2 + (2*c1*kon)/koff + (c2*kon)/koff2))/kf2^2 + 
               (koff2*(kf2*(1 + (c1*kon)/koff)^2 - (c1*kf*kon*(3 + (3*c1*kon)/
                      koff + (c2*kon)/koff2))/koff))/kf2))/kf + 
            (koff*koff2*((c1*c2*(1 + koff2/kf2)*kon^2)/koff + kf*
                (1 + (c1*kon)/koff + (c1^2*kon^2)/koff^2 - (c1*c2*kon^2)/
                  (koff*koff2) + (koff2*(2 + (c2*kon)/koff2 + (2*c1^2*kon^2)/
                     koff^2 - (c1*kon*(-2 + (c2*kon)/koff2))/koff))/kf2 + 
                 (koff2^2*((c1^2*kon^2)/koff^2 + (c1*kon*(1 + (c2*kon)/
                        koff2))/koff + (1 + (c2*kon)/koff2)^2))/kf2^2)))/kf))/
          (koff*koff2^2) - (2*c1*(1 + koff2/kf2)^2*kon*
           ((c2*koff*(1 + koff/kf)*kon*(1 + koff/kf + (c1*kon)/kf))/koff2 + 
            (c1*koff2^2*kon*(1 + (c2*kon)/koff2 + (koff*(2 + (c1*kon)/koff + 
                  (2*c2*kon)/koff2))/kf))/(kf2*koff) + 
            (koff2*(-((c1*c2*(1 + koff/kf)*kon^2)/koff2) + kf2*
                ((c1^2*kon^2)/(kf*koff) + (c2*koff*(1 + koff/kf)*kon*
                   (1 + (c2*kon)/koff2))/(kf*koff2) + (c1*kon*(1 + 
                    (koff*(2 + (c2*kon)/koff2))/kf))/koff)))/kf2))/
          (koff*koff2) + (2*c1*(1 + koff2/kf2)^2*kon*
           (-(c2*koff2*kon*(1 + (c2*kon)/koff2)) + 
            (koff*koff2*(-(c2*kon*(3 + (c1*kon)/koff + (3*c2*kon)/koff2)) + 
               kf*(1 + (c2*kon)/koff2 + (c2^2*kon^2)/koff2^2 - (c1*c2*kon^2)/
                  (koff*koff2) + (koff2*(1 + (c2*kon)/koff2)^2)/kf2)))/kf + 
            (koff^3*((c1*c2*kf*kon^2)/(koff*koff2) + (koff2^2*
                 (1 + (c1*kon)/koff + (c2*kon)/koff2)^2)/kf2 + (koff2*
                 ((c1*c2*kf*kon^2)/(koff*koff2) + kf2*(1 + (c2*kon)/koff2 + 
                    (c1^2*kon^2)/koff^2 + (c2^2*kon^2)/koff2^2 + 
                    (c1*kon*(2 + (c2*kon)/koff2))/koff)))/kf2))/kf^2 + 
            (koff^2*((c1*c2*kf^2*kon^2)/(koff*koff2) + (koff2^2*
                 (2 + (c1*kon)/koff + (2*c2*kon)/koff2)*(kf + (c2*kf*kon)/
                   koff2 - (c2*kf2*kon)/koff2))/kf2 + (kf*koff2*
                 ((c1*c2*kf*kon^2)/(koff*koff2) + kf2*((c1*kon)/koff - 
                    (c1*c2*kon^2)/(koff*koff2) + 2*(1 + (c2*kon)/koff2 + 
                      (c2^2*kon^2)/koff2^2))))/kf2))/kf^2))/(koff^2*koff2) + 
         c2*(1 + koff/kf)^2*(1 + koff2/kf2)*kon*(1 + (c1*kon)/koff + 
           (c2*kon)/koff2)*((c1*kon)/koff + (c2*kon)/koff2 + 
           (c2*koff*kon)/(kf*koff2) + (c1*koff2*kon)/(kf2*koff))*t + 
         c1*(1 + koff/kf)*(1 + koff2/kf2)^2*kon*(1 + (c1*kon)/koff + 
           (c2*kon)/koff2)*((c1*kon)/koff + (c2*kon)/koff2 + 
           (c2*koff*kon)/(kf*koff2) + (c1*koff2*kon)/(kf2*koff))*t + 
         (c2*(1 + koff/kf)^2*(1 + koff2/kf2)*kon*(1 + (c1*kon)/koff + 
            (c2*kon)/koff2)*(c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*
                 kon)/koff))/kf)*t)/koff2 + (c1*(1 + koff/kf)*
           (1 + koff2/kf2)^2*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)*
           (c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*kon)/koff))/kf)*
           t)/koff - 2*(1 + koff/kf)*(1 + koff2/kf2)*(1 + (c1*kon)/koff + 
           (c2*kon)/koff2)*((c1*kon)/koff + (c2*kon)/koff2 + 
           (c2*koff*kon)/(kf*koff2) + (c1*koff2*kon)/(kf2*koff))*
          (c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*kon)/koff))/kf)*
          t))/(2*(1 + koff/kf)^3*(1 + koff2/kf2)^3*
        (1 + (c1*kon)/koff + (c2*kon)/koff2)^3), 
      (t*((c1*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)^2)/(1 + koff/kf) + 
         (c2*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)^2)/(1 + koff2/kf2) - 
         (2*c1*kon*(-((c1*c2*koff2*kon^2)/koff) + 
            (c2*koff^3*kon*(-((c1*kf*kon)/koff) + (koff2*(kf2 - (c1*kf*kon)/
                   koff + (c1*kf2*kon)/koff))/kf2))/(kf^2*koff2) + 
            (koff*koff2*(c2*kon*(1 - (c1*kon)/koff + (c2*kon)/koff2) + kf*
                ((c2*kon)/koff2 + (c1*kon*(1 + koff2/kf2 + (c2*kon)/kf2 + 
                    (2*c2*kon)/koff2))/koff)))/kf + 
            (koff^2*(-((c1*c2*kf^2*kon^2)/(koff*koff2)) + (koff2^2*
                 ((c2*kf2*kon*(1 + (c2*kon)/koff2))/koff2 + (c1*kf*kon*
                    (2 + (c1*kon)/koff + (2*c2*kon)/koff2))/koff))/kf2 + 
               (kf*koff2*(-((c1*c2*kf*kon^2)/(koff*koff2)) + 
                  kf2*((2*c1*kon)/koff + (2*c2*kon)/koff2 + (c1^2*kon^2)/
                     koff^2 + (4*c1*c2*kon^2)/(koff*koff2))))/kf2))/kf^2))/
          (koff*(1 + koff/kf)^3*koff2*(1 + koff2/kf2)) - 
         (2*c2*kon*(-((c1*c2*koff2*(1 + koff2/kf2)*kon^2)/koff) + 
            (koff^2*(-((c1*c2*kf*kon^2)/(koff*koff2)) + (koff2*
                 ((c2*kf2*kon*(1 + (c1*kon)/koff))/koff2 + (c1*kf*kon*
                    (1 + (c1*kon)/koff - (c2*kon)/koff2))/koff))/kf2 + 
               (koff2^2*((c1*kf*kon*(1 + (c1*kon)/koff))/koff + 
                  (c2*kf2*kon*(2 + (2*c1*kon)/koff + (c2*kon)/koff2))/koff2))/
                kf2^2))/kf + (koff*koff2*(-((c1*c2*(1 + koff2/kf2)*kon^2)/
                 koff) + kf*((c2*kon*(1 + (koff2*(2 + (c2*kon)/koff2))/kf2))/
                  koff2 + (c1*kon*(1 + (2*c2*kon)/koff2 + (koff2^2*
                      (1 + (c2*kon)/koff2))/kf2^2 + (koff2*(2 + (4*c2*kon)/
                        koff2))/kf2))/koff)))/kf))/(koff*(1 + koff/kf)*koff2*
           (1 + koff2/kf2)^3) + (c2*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)*
           (c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*kon)/koff))/kf)*
           t)/((1 + koff/kf)*(1 + koff2/kf2)^2) + 
         (c1*kon*(1 + (c1*kon)/koff + (c2*kon)/koff2)*
           (c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*kon)/koff))/kf)*
           t)/((1 + koff/kf)^2*(1 + koff2/kf2)) - 
         ((1 + (c1*kon)/koff + (c2*kon)/koff2)*
           (c2*kon + (koff*(c2*kon + (c1*kf*(1 + koff2/kf2)*kon)/koff))/kf)^2*
           t)/((1 + koff/kf)^2*(1 + koff2/kf2)^2)))/
       (1 + (c1*kon)/koff + (c2*kon)/koff2)^3}}
