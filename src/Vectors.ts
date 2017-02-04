namespace Vectors {
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
        copy(v: Vector) {
            this.x = v.x;
            this.y = v.y;
            this.z = v.z;
            return this;
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

        static oadd(u: Vector, v: Vector, out: Vector): Vector {
            return out.set(u.x+v.x, u.y+v.y, u.z+v.z);
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
            return this;
        }

        isub(v: Vector) {
            this.x -= v.x;
            this.y -= v.y;
            this.z -= v.z;
            return this;
        }

        imult(s: number) {
            this.x *= s;
            this.y *= s;
            this.z *= s;
            return this;
        }

        izero() {
            this.x = this.y = this.z = 0;
            return this;
        }

        free() {
            VectorPool.uptake(this);
        }

        static getVector(): Vector;
        static getVector(n: number): Vector[];
        static getVector(n?: number): any {
            if (n) {
                var vs: Vector[] = [];
                while ((n--) > 0)
                    vs.push(new Vector());
                return vs;
            }
            return (new Vector());
        }

    }

    export class VectorPool {
        public static counter: number = 0;
        private static allocated: Vector[] = [];
        private static requested: Vector[] = [];
        
        constructor() {}

        static collect() {
            var v: Vector;
            while (VectorPool.requested.length > 0) {
                v = VectorPool.requested.pop();
                VectorPool.allocated.push(v);
            }
        }
        static uptake(v: Vector) {
            this.allocated.push(v);
        }


        static getVector(): Vector;
        static getVector(n: number): Vector[];
        static getVector(n?: number): any {
            if (n) {
                var vs: Vector[] = [];
                while ((n--) > 0)
                    vs.push(VectorPool.getVector());
                return vs;
            }
            else {
                if (VectorPool.allocated.length===0) {
                    VectorPool.allocated.push(new Vector());
                    VectorPool.counter++;
                    // console.log("allocating new vector");
                    
                }
                var v = VectorPool.allocated.pop()
                VectorPool.requested.push(v);
                v.izero();
                return v;
            }
        }

    }

    //---------------------------------------------------
    // export class VPool {
    //     public static counter: number = 0;
    //     private allocated: Vector[] = [];
    //     private requested: Vector[] = [];
        
    //     constructor() {}

    //     collect() {
    //         var v: Vector;
    //         while (this.requested.length > 0) {
    //             v = this.requested.pop();
    //             this.allocated.push(v);
    //         }
    //     }

    //     getVector(): Vector;
    //     getVector(n: number): Vector[];
    //     getVector(n?: number): any {
    //         if (n) {
    //             var vs: Vector[] = [];
    //             while ((n--) > 0)
    //                 vs.push(VectorPool.getVector());
    //             return vs;
    //         }
    //         else {
    //             if (this.allocated.length===0) {
    //                 this.allocated.push(new Vector());
    //                 VPool.counter++;
    //             }
    //             var v = this.allocated.pop()
    //             this.requested.push(v);
    //             v.izero();
    //             return v;
    //         }
    //     }

    //     free() {
    //         Pool.uptake(this);
    //     }
    // }


    // export class Pool {
    //     private static counter: number = 0;
    //     private static allocated: VPool[] = [];
    //     private static requested: VPool[] = [];
        
    //     constructor() {}

    //     static uptake(v: VPool) {
    //         var i = Pool.requested.indexOf(v);
    //         Pool.requested.splice(i);
    //         Pool.allocated.push(v);
    //     }

    //     static getVPool(): VPool {
    //         if (Pool.allocated.length===0) {
    //             Pool.allocated.push(new VPool());
    //             Pool.counter++;
    //         }
    //         var v = Pool.allocated.pop();
    //         Pool.requested.push(v);
    //         return v;
    //     }
    // }
}