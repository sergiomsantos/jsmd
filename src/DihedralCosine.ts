/// <reference path="ForceField.ts" />

namespace ForceField {
    export class DihedralCosine implements IForceFieldComponent {
        a1: number;
        a2: number;
        a3: number;
        a4: number;
        kf: number;
        pn: number;
        ph: number;
        
        constructor(a1: number, a2: number, a3: number, a4: number, kf: number, pn: number, ph: number) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
            this.a4 = a4;
            this.kf = kf;
            this.pn = pn;
            this.ph = ph;
        }

        static eval(dihedrals: DihedralCosine[], coords: Vector[], forces?: Vector[]): number {
            var v12: Vector = VectorPool.getVector();
            var v32: Vector = VectorPool.getVector();
            var v34: Vector = VectorPool.getVector();
            var tmp: Vector = VectorPool.getVector();
            var M: Vector = VectorPool.getVector();
            var N: Vector = VectorPool.getVector();
            var f1: Vector = VectorPool.getVector();
            var f2: Vector = VectorPool.getVector();
            var f3: Vector = VectorPool.getVector();
            var f4: Vector = VectorPool.getVector();
            var d32Sq: number,
                d32: number,
                phi: number,
                dVdphi_x_d32: number,
                energy: number = 0.0;

            dihedrals.forEach(dihd => {
                Vector.osub(coords[dihd.a1], coords[dihd.a2], v12);
                Vector.osub(coords[dihd.a3], coords[dihd.a2], v32);
                Vector.osub(coords[dihd.a3], coords[dihd.a4], v34);
                Vector.ocross(v12,v32,M);
                Vector.ocross(v32,v34,N);
                d32Sq = v32.normSq();
                d32 = Math.sqrt(d32Sq);
                phi = Math.atan2(d32*v12.dot(N), M.dot(N)) - Math.PI;
                
                if (forces) {
                    dVdphi_x_d32 = dihd.kf * dihd.pn * Math.sin(dihd.pn*phi) * d32;
                    Vector.omult(M, -dVdphi_x_d32/M.dot(M), f1);    // F1
                    Vector.omult(N,  dVdphi_x_d32/N.dot(N), f4);    // F4
                    
                    f3.set(-f4.x, -f4.y, -f4.z);
                    f3.isub(Vector.omult(f1, v12.dot(v32)/d32Sq, tmp));
                    f3.iadd(Vector.omult(f4, v34.dot(v32)/d32Sq, tmp));
                    
                    tmp.set(-f1.x-f3.x-f4.x, -f1.y-f3.y-f4.y, -f1.z-f3.z-f4.z);
                    forces[dihd.a1].iadd(f1);
                    forces[dihd.a2].iadd(tmp);
                    forces[dihd.a3].iadd(f3);
                    forces[dihd.a4].iadd(f4);
                }
                energy += dihd.kf * (1.0 + Math.cos(dihd.pn*phi - dihd.ph));
            });

            v12.free();
            v32.free();
            v34.free();
            tmp.free();
            M.free();
            N.free();
            f1.free();
            f2.free();
            f3.free();
            f4.free();
            // VectorPool.collect();
            return energy;
        }
    }
}