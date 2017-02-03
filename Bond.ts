/// <reference path="ForceField.ts" />

namespace ForceField {

    export class Bond implements IForceFieldComponent {
        a1: number;
        a2: number;
        kf: number;
        r0: number;
        constructor(a1: number, a2: number, kf: number, r0: number) {
            this.a1 = a1;
            this.a2 = a2;
            this.kf = kf;
            this.r0 = r0;
        }
        static eval(bonds: Bond[], coords: Vector[], forces?: Vector[]): number {
            var v12: Vector = VectorPool.getVector();// new Vector();
            var d12: number,
                dr: number,
                energy: number = 0.0;
            
            bonds.forEach(bond => {
                Vector.osub(coords[bond.a2], coords[bond.a1], v12);
                d12 = v12.norm();
                dr = d12-bond.r0;
                if (forces) {
                    v12.imult(bond.kf*dr/d12);
                    forces[bond.a1].iadd(v12);
                    forces[bond.a2].isub(v12);
                }
                energy += bond.kf * dr * dr;
            });

            v12.free();
            // VectorPool.collect();
            return 0.5*energy;
        }

    }
}