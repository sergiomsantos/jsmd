"use strict";
class VectorPool {
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



export class Vector {
    x: number;
    y: number;
    z: number;
    // private static vlist: Vector[] = [];
    // private static vrequested: Vector[] = [];
    
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
    
    // static getVector(): Vector {
    //     if (Vector.vlist.length===0) {
    //         Vector.vlist.push(new Vector());
    //     }
    //     var v = Vector.vlist.pop()
    //     Vector.vrequested.push(v);
    //     v.set(0,0,0);
    //     return v;
    // }
    // // static putVector(v: Vector) {
    // //     Vector.vlist. push(v);
    // // }

    // static collect() {
    //     var v: Vector;
    //     while (Vector.vrequested.length > 0) {
    //         v = Vector.vrequested.pop();
    //         Vector.vlist.push(v);
    //     }
    // }

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
    // toString(): string {
    //     return `v={${this.x},${this.y},${this.z}}`
    // }
    // release() {
    //     Vector.vlist.push(this);
    // }

}

interface IForceFieldBondedComponent {
    eval(coords: Vector[], forces?: Vector[]): number;
}

export class Bond implements IForceFieldBondedComponent {
    a1: number;
    a2: number;
    kf: number;
    r0: number;
    constructor(_a1: number, _a2: number, _kf: number, _r0: number) {
        this.a1 = _a1;
        this.a2 = _a2;
        this.kf = _kf;
        this.r0 = _r0;
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
        VectorPool.collect();
        //Vector.putVector(v12);
        return 0.5*energy;
    }

    eval(coords: Vector[], forces?: Vector[]): number {
        return 0;
        // var v12 = Vector.sub(coords[this.a2], coords[this.a1]);
        // var d12 = v12.norm();
        // var dr = d12-this.r0;
        // if (forces) {
        //     v12.imult(this.kf*dr/d12);
        //     forces[this.a1].iadd(v12);
        //     forces[this.a2].isub(v12);
        // }
        // return this.kf * dr * dr;
    }
}

export class Angle implements IForceFieldBondedComponent {
    a1: number;
    a2: number;
    a3: number;
    kf: number;
    r0: number;
    constructor(_a1: number, _a2: number, _a3: number, _kf: number, _r0: number) {
        this.a1 = _a1;
        this.a2 = _a2;
        this.a3 = _a3;
        this.kf = _kf;
        this.r0 = _r0;
    }

    eval(coords: Vector[], forces?: Vector[]) {
        var v12 = Vector.sub(coords[this.a2], coords[this.a1]);
        var v32 = Vector.sub(coords[this.a2], coords[this.a3]);
        var d12Sq = v12.normSq();
        var d32Sq = v32.normSq();
        var d12xd32 = Math.sqrt(d12Sq*d32Sq);

        var ctheta = v12.dot(v32) / d12xd32;
        var dtheta = Math.acos(ctheta) - this.r0;

        if (forces) {
            var fc = this.kf*dtheta / Math.sqrt(1.0 - ctheta);
            var f1 = new Vector(fc*(ctheta*v12.x/d12Sq - v32.x/d12xd32),
                                fc*(ctheta*v12.y/d12Sq - v32.y/d12xd32),
                                fc*(ctheta*v12.z/d12Sq - v32.z/d12xd32));
            var f3 = new Vector(fc*(ctheta*v32.x/d32Sq - v12.x/d12xd32),
                                fc*(ctheta*v32.y/d32Sq - v12.y/d12xd32),
                                fc*(ctheta*v32.z/d32Sq - v12.z/d12xd32));
            forces[this.a1].iadd(f1);
            forces[this.a3].iadd(f3);
            f1.iadd(f3);
            forces[this.a2].isub(f1);
        }

        return 0.5 * this.kf * dtheta * dtheta;
    }
}




class DihedralRB implements IForceFieldBondedComponent {
    a1: number;
    a2: number;
    a3: number;
    a4: number;
    c0: number;
    c1: number;
    c2: number;
    c3: number;
    c4: number;
    c5: number;
    constructor(_a1: number, _a2: number, _a3: number, _a4: number,
                _c0: number, _c1: number, _c2: number, _c3: number, _c4: number, _c5: number) {
        this.a1 = _a1;
        this.a2 = _a2;
        this.a3 = _a3;
        this.a4 = _a4;
        this.c0 = _c0;
        this.c1 = _c1;
        this.c2 = _c2;
        this.c3 = _c3;
        this.c4 = _c4;
        this.c5 = _c5;
    }

    static eval(dihedrals: DihedralRB[], coords: Vector[], forces?: Vector[]): number {
        var v12: Vector = VectorPool.getVector();
        var v32: Vector = VectorPool.getVector();
        var v34: Vector = VectorPool.getVector();
        var tmp: Vector = VectorPool.getVector();
        var  f1: Vector = VectorPool.getVector();
        var  f2: Vector = VectorPool.getVector();
        var  f3: Vector = VectorPool.getVector();
        var  f4: Vector = VectorPool.getVector();
        var   M: Vector = VectorPool.getVector();
        var   N: Vector = VectorPool.getVector();
        var d32Sq: number,
            d32: number,
            phi: number,
            cpsi: number,
            cphi: number,
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
            phi = Math.atan2(d32*v12.dot(N), M.dot(N)) ;
            cpsi = Math.cos(phi-Math.PI);
            energy += dihd.c0 + cpsi*(dihd.c1 + cpsi*(dihd.c2 + cpsi*(dihd.c3 + cpsi*(dihd.c4 + cpsi*dihd.c5))));

            if (forces) {
                cphi = -cpsi; // cos(phi-PI) = cos(psi) = -cos(phi)
                dVdphi_x_d32 = Math.sin(phi) * d32 * (dihd.c1 + cphi*(-2.0*dihd.c2 + cphi*(3.0*dihd.c3 + cphi*(-4.0*dihd.c4 + cphi*5.0*dihd.c5))));

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
        });
        VectorPool.collect();
        // v12.release();
        // v32.release();
        // v34.release();
        // tmp.release();
        // f1.release();
        // f2.release();
        // f3.release();
        // f4.release();
        // M.release();
        // N.release();
        return energy;
    }

    eval(coords: Vector[], forces?: Vector[]): number {
        return 0;
    }
    // eval(coords: Vector[], forces?: Vector[]) {
    //     var v12 = Vector.sub(coords[this.a1], coords[this.a2]);
    //     var v32 = Vector.sub(coords[this.a3], coords[this.a2]);
    //     var v34 = Vector.sub(coords[this.a3], coords[this.a4]);
    //     var M = Vector.cross(v12,v32);
    //     var N = Vector.cross(v32,v34);
    //     var d32Sq = v32.normSq();
    //     var d32 = Math.sqrt(d32Sq);

    //     var phi = Math.atan2(d32*v12.dot(N), M.dot(N)) - Math.PI;

    //     var d12 = v12.norm();
    //     var dr = d12-this.r0;
    //     if (forces) {
    //         var dVdphi_x_d32 = this.kf * this.pn * Math.sin(this.pn*phi) * d32;
    //         M.imult(-dVdphi_x_d32/M.dot(M));    // F1
    //         N.imult( dVdphi_x_d32/N.dot(N));    // F4
    //         var f3 = new Vector(-N.x, -N.y, -N.z);
    //         f3
    //         forces[this.a1].iadd(M);
    //         forces[this.a3].iadd(f3);
    //         forces[this.a4].iadd(N);

    //         forces[this.a2].iadd(f1);
    //     }
    //     return this.kf * (1.0 + Math.cos(this.pn*phi - this.ph));
    // }
}



class Dihedral implements IForceFieldBondedComponent {
    a1: number;
    a2: number;
    a3: number;
    a4: number;
    kf: number;
    pn: number;
    ph: number;
    constructor(_a1: number, _a2: number, _a3: number, _a4: number, _kf: number, _pn: number, _ph: number) {
        this.a1 = _a1;
        this.a2 = _a2;
        this.a3 = _a3;
        this.a4 = _a4;
        this.kf = _kf;
        this.pn = _pn;
        this.ph = _ph;
    }

    static eval(dihedrals: Dihedral[], coords: Vector[], forces?: Vector[]): number {
        var v12: Vector = VectorPool.getVector();//new Vector();
        var v32: Vector = VectorPool.getVector();//new Vector();
        var v34: Vector = VectorPool.getVector();//new Vector();
        var tmp: Vector = VectorPool.getVector();//new Vector();
        var M: Vector = VectorPool.getVector();//new Vector();
        var N: Vector = VectorPool.getVector();//new Vector();
        var f1: Vector = VectorPool.getVector();//new Vector();
        var f2: Vector = VectorPool.getVector();//new Vector();
        var f3: Vector = VectorPool.getVector();//new Vector();
        var f4: Vector = VectorPool.getVector();//new Vector();
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
            // console.log('phi', phi);
            
            if (forces) {
                dVdphi_x_d32 = dihd.kf * dihd.pn * Math.sin(dihd.pn*phi) * d32;
                Vector.omult(M, -dVdphi_x_d32/M.dot(M), f1);    // F1
                // console.log(f1,M,dVdphi_x_d32);
                
                Vector.omult(N,  dVdphi_x_d32/N.dot(N), f4);    // F4
                f3.set(-f4.x, -f4.y, -f4.z);
                f3.isub(Vector.omult(f1, v12.dot(v32)/d32Sq, tmp));
                f3.iadd(Vector.omult(f4, v34.dot(v32)/d32Sq, tmp));
                // console.log(f1, f4);
                
                tmp.set(-f1.x-f3.x-f4.x, -f1.y-f3.y-f4.y, -f1.z-f3.z-f4.z);
                forces[dihd.a1].iadd(f1);
                forces[dihd.a2].iadd(tmp);
                forces[dihd.a3].iadd(f3);
                forces[dihd.a4].iadd(f4);
            }
            energy += dihd.kf * (1.0 + Math.cos(dihd.pn*phi - dihd.ph));
        });
        VectorPool.collect();
        // Vector.putVector(v12);
        // Vector.putVector(v32);
        // Vector.putVector(v34);
        // Vector.putVector(tmp);
        // Vector.putVector(f1);
        // Vector.putVector(f2);
        // Vector.putVector(f3);
        // Vector.putVector(f4);
        // Vector.putVector(M);
        // Vector.putVector(N);
        return energy;
    }

    eval(coords: Vector[], forces?: Vector[]): number {
        return 0;
    }
    // eval(coords: Vector[], forces?: Vector[]) {
    //     var v12 = Vector.sub(coords[this.a1], coords[this.a2]);
    //     var v32 = Vector.sub(coords[this.a3], coords[this.a2]);
    //     var v34 = Vector.sub(coords[this.a3], coords[this.a4]);
    //     var M = Vector.cross(v12,v32);
    //     var N = Vector.cross(v32,v34);
    //     var d32Sq = v32.normSq();
    //     var d32 = Math.sqrt(d32Sq);

    //     var phi = Math.atan2(d32*v12.dot(N), M.dot(N)) - Math.PI;

    //     var d12 = v12.norm();
    //     var dr = d12-this.r0;
    //     if (forces) {
    //         var dVdphi_x_d32 = this.kf * this.pn * Math.sin(this.pn*phi) * d32;
    //         M.imult(-dVdphi_x_d32/M.dot(M));    // F1
    //         N.imult( dVdphi_x_d32/N.dot(N));    // F4
    //         var f3 = new Vector(-N.x, -N.y, -N.z);
    //         f3
    //         forces[this.a1].iadd(M);
    //         forces[this.a3].iadd(f3);
    //         forces[this.a4].iadd(N);

    //         forces[this.a2].iadd(f1);
    //     }
    //     return this.kf * (1.0 + Math.cos(this.pn*phi - this.ph));
    // }
}



class Atom {
    id: number;
    name: string;
    sigma: number;
    charge: number;
    epsilon: number;
    blist: number[];
    constructor(id: number, name: string, sigma: number, epsilon: number, charge: number) {
        this.sigma = sigma;
        this.epsilon = epsilon;
        this.charge = charge;
        this.id = id;
        this.name = name;
        this.blist = [Number.POSITIVE_INFINITY];
    }
    addBondedAtom(atomId: number) {
        if (this.blist.indexOf(atomId) < 0) {
            this.blist.push(atomId);
        }
        this.blist.sort((a,b) => {return a-b;});

    }
    static eval(atoms: Atom[], coords: Vector[], forces?: Vector[]): number {
        var energy: number = 0.0;
        var nAtoms = atoms.length;
        var bidx, bid;
        var atomi: Atom, atomj: Atom;
        var vij: Vector = VectorPool.getVector();
        var dijSq: number, sij: number, eij: number;
        var fc: number, lj6: number;
        var Eq: number, Elj: number;

        for (var i=0; i<nAtoms-1; i++) {
            atomi = atoms[i];
            bidx = 0;
            bid = atomi[bidx];
            
            
            for (var j=i+1; j<nAtoms; j++) {
                if (j===bid) {
                    bidx++;
                    bid = atomi[bidx];
                    continue;
                }
                Vector.osub(coords[j], coords[i], vij);
                dijSq = vij.normSq();

                atomj = atoms[j];
                //console.log(i,atomi,j,atomj);
                sij = atomi.sigma + atomj.sigma;
                eij = atomi.epsilon * atomj.epsilon;
                lj6 = Math.pow(sij*sij/dijSq, 3);
                
                Elj = 4.0 * eij * (lj6*lj6-lj6);
                Eq = atomi.charge * atomj.charge / Math.sqrt(dijSq);
                //console.log(Elj, Eq);
                
                if (forces) {
                    fc = -(24*eij*(2*lj6*lj6-lj6) + Eq)/dijSq;
                    
                    vij.imult(fc);
                    //console.log(vij);
                    forces[i].iadd(vij);
                    forces[j].isub(vij);
                }

                energy += Elj + Eq;
            }
        }
        VectorPool.collect();
        //Vector.putVector(vij);
        return energy;
    }
}




class System {
    private bonded: IForceFieldBondedComponent[] = [];
    private nonbonded: IForceFieldBondedComponent[] = [];
    private dihedrals: Dihedral[] = [];
    private angles: Angle[] = [];
    private bonds: Bond[] = [];
    private atoms: Atom[] = [];

    constructor() {}
    newBond(a1: number, a2: number, kf: number, r0: number): Bond {
        var item = new Bond(a1, a2, kf, r0);
        //this.bonded.push(item);
        this.bonds.push(item);
        return item;
    }

    newAngle(a1: number, a2: number, a3: number, kf: number, r0: number): Bond {
        var item = new Angle(a1, a2, a3, kf, r0);
        this.bonded.push(item);
        this.angles.push(item);
        return item;
    }

    newDihedral(a1: number, a2: number, a3: number, a4: number, 
                kf: number, pn: number, ph: number): Dihedral {
        var item = new Dihedral(a1, a2, a3, a4, kf, pn, ph);
        //this.bonded.push(item);
        this.dihedrals.push(item);
        return item;
    }
    newAtom(id: number, name: string, sigma: number, epsilon: number, charge: number): Atom {
        var item = new Atom(id, name, sigma, epsilon, charge);
        this.atoms.push(item);
        return item;
    }

    eval(coords: Vector[], forces?: Vector[]): number {
        var energy = 0;
        // this.bonded.forEach(element => {
        //     energy += element.eval(coords, forces)
        // });
        // energy +=     Bond.eval(this.bonds,     coords, forces);
        // energy += Dihedral.eval(this.dihedrals, coords, forces);
        energy +=     Atom.eval(this.atoms,     coords, forces);
        return energy;
    }

}

var d=1.5;
var coords = [
    new Vector(0,0,0),
    new Vector(d,0,0),
    new Vector(d,d,0),
    new Vector(d,d,d)
];

var forces = [
    new Vector(),
    new Vector(),
    new Vector(),
    new Vector()
];


// for (let v of coords)
//     console.log(v);

// for (let v of forces)
//     console.log(v);

// var bond = new Bond(0, 1, 10, 4);
// console.log(bond);

// console.log(bond.eval(coords));
// for (let v of forces)
//     console.log(v);
// console.log(bond.eval(coords, forces));
// for (let v of forces)
//     console.log(v);

// var angle = new Angle(0, 1, 2, 10, Math.PI);
// console.log(angle);

// console.log(angle.eval(coords));
// for (let v of forces)
//     console.log(v);
// console.log(angle.eval(coords, forces));
// for (let v of forces)
//     console.log(v);


var top = {
    bonds: [
        {
            a1: 0,
            a2: 1,
            kf: 10,
            r0: 5
        },
        {
            a1: 1,
            a2: 2,
            kf: 10,
            r0: 5
        },
        {
            a1: 2,
            a2: 3,
            kf: 10,
            r0: 5
        }
    ],
    angles: [
        {
            a1: 0,
            a2: 1,
            a3: 2,
            kf: 10,
            r0: 2*Math.PI/3
        },
        {
            a1: 1,
            a2: 2,
            a3: 3,
            kf: 10,
            r0: 2*Math.PI/3
        }
    ],
    dihedrals: [
        {a1: 0, a2: 1, a3: 2, a4: 3, kf: 10, pn:3, ph:2*Math.PI/3}
    ],
    atoms: [
        {id:1, name:'C1', sigma: 1, epsilon: 0, charge: -1},
        {id:2, name:'C2', sigma: 1, epsilon: 0, charge:  1},
        {id:3, name:'C3', sigma: 1, epsilon: 0, charge: -1},
        {id:4, name:'C4', sigma: 1, epsilon: 0, charge:  1},
    ]
}

var system = new System();
top.bonds.forEach(item => {
    system.newBond(item.a1, item.a2, item.kf, item.r0);
});

top.angles.forEach(item => {
    system.newAngle(item.a1, item.a2, item.a3, item.kf, item.r0);
});

top.dihedrals.forEach(item => {
    system.newDihedral(item.a1, item.a2, item.a3, item.a4, item.kf, item.pn, item.ph);
});

top.atoms.forEach(item => {
    var atom = system.newAtom(item.id, item.name, item.sigma, item.epsilon, item.charge);
    // console.log(atom);
});


// console.log(system.eval(coords));
// for (let v of forces)
//     console.log(v);

// console.log(system.eval(coords, forces));
// for (let v of forces)
//     console.log(v);

console.log('   ', 4);
console.log('frame', e_new);
for (let v of coords)
    console.log('C', v.x, v.y, v.z);

var e_old = 10;
var e_new = 0;
while (Math.abs(e_new - e_old) > 0.0001) {
    forces.forEach(element => {
        element.izero();
    });

    e_old = e_new;
    e_new = system.eval(coords, forces);
    // console.log(e_new);
    
    for (var i=0; i<coords.length; i++) {
        forces[i].imult(0.001);
        coords[i].iadd(forces[i]);
    }
    console.log('   ', 4);
    console.log('frame', e_new);
    for (let v of coords)
    console.log('C', v.x, v.y, v.z);
    
}

// var e = system.eval(coords, forces);
    // console.log(e);
// for (let v of forces)
//     console.log(v);
// for (let v of coords)
//     console.log(v);