/// <reference path="ForceField.ts" />

namespace ForceField {

    export class Angle implements IForceFieldComponent {
    
        a1: number;
        a2: number;
        a3: number;
        kf: number;
        r0: number;
        constructor(a1: number, a2: number, a3: number, kf: number, r0: number) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
            this.kf = kf;
            this.r0 = r0;
        }

        static eval(angles: Angle[], coords: Vector[], forces?: Vector[]) {
            var v12: Vector = VectorPool.getVector();
            var v32: Vector = VectorPool.getVector();
            var f1: Vector = VectorPool.getVector();
            var f3: Vector = VectorPool.getVector();
            var d32Sq: number,
                d12Sq: number,
                d12xd32: number,
                fc: number,
                energy: number = 0.0;

            angles.forEach(angle => {
                Vector.osub(coords[angle.a2], coords[angle.a1], v12);
                Vector.osub(coords[angle.a2], coords[angle.a3], v32);
                d12Sq = v12.normSq();
                d32Sq = v32.normSq();
                d12xd32 = Math.sqrt(d12Sq*d32Sq);
                
                var ctheta = v12.dot(v32) / d12xd32;
                var dtheta = Math.acos(ctheta) - angle.r0;

                if (forces) {
                    fc = angle.kf*dtheta / Math.sqrt(1.0 - ctheta);
                    f1.set(fc*(ctheta*v12.x/d12Sq - v32.x/d12xd32),
                           fc*(ctheta*v12.y/d12Sq - v32.y/d12xd32),
                           fc*(ctheta*v12.z/d12Sq - v32.z/d12xd32));
                    f3.set(fc*(ctheta*v32.x/d32Sq - v12.x/d12xd32),
                           fc*(ctheta*v32.y/d32Sq - v12.y/d12xd32),
                           fc*(ctheta*v32.z/d32Sq - v12.z/d12xd32));
                    forces[angle.a1].iadd(f1);
                    forces[angle.a3].iadd(f3);
                    f1.iadd(f3);
                    forces[angle.a2].isub(f1);
                }
                energy += angle.kf * dtheta * dtheta;
            });
            
            VectorPool.collect();
            return 0.5*energy;
        }
    }
}