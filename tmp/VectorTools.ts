namespace VectorTools {
    export class Vector {
        x: number;
        y: number;
        z: number;

        constructor(_x: number=0, _y:number=0, _z:number=0) {
            this.x = _x;
            this.y = _y;
            this.z = _z;        
        }

        dot(other: Vector): number {
            return (this.x*other.x + this.y*other.y + this.z*other.z);
        }

        normSq(): number {
            return this.dot(this);
        }

        norm(): number {
            return Math.sqrt(this.normSq())
        }

        set(x: number, y: number, z: number) {
            this.x = x;
            this.y = y;
            this.z = z;
            return this;
        }

        static sub(v1: Vector, v2: Vector): Vector {
            return new Vector(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z)
        }

        static osub(u: Vector, v: Vector, out: Vector): Vector {
            return out.set(u.x-v.x, u.y-v.y, u.z-v.z);
        }

        static add(u: Vector, v: Vector): Vector {
            return new Vector(u.x+v.x, u.y+v.y, u.z+v.z)
        }

        static cross(u: Vector, v: Vector): Vector {
            return new Vector(u.y*v.z - u.z*v.y,
                            u.z*v.x - u.x*v.z,
                            u.x*v.y - u.y*v.x);
        }

        static ocross(u: Vector, v: Vector, out: Vector): Vector {
            return out.set(u.y*v.z - u.z*v.y,
                        u.z*v.x - u.x*v.z,
                        u.x*v.y - u.y*v.x);
        }

        static mult(v: Vector, s: number): Vector {
            return new Vector(v.x*s, v.y*s, v.z*s);
        }

        static omult(v: Vector, s: number, out: Vector) {
            return out.set(s*v.x, s*v.y, s*v.z);
        }

        iadd(v: Vector) {
            this.x += v.x;
            this.y += v.y;
            this.z += v.z;
        }

        isub(v: Vector) {
            this.x -= v.x;
            this.y -= v.y;
            this.z -= v.z;
        }

        imult(s: number) {
            this.x *= s;
            this.y *= s;
            this.z *= s;
        
        }

        izero() {
            this.x = this.y = this.z = 0;
        }
    }

    export class VectorPool {
        private static available: Vector[] = [];
        private static requested: Vector[] = [];
        constructor() {}

        static collect() {
            var v: Vector;
            while (VectorPool.requested.length > 0) {
                v = VectorPool.requested.pop();
                VectorPool.available.push(v);
            }
        }

        static getVector(): Vector {
            if (VectorPool.available.length===0) {
                VectorPool.available.push(new Vector());
            }
            var v = VectorPool.available.pop()
            VectorPool.requested.push(v);
            v.set(0,0,0);
            return v;
        }

}






}