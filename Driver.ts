/// <reference path="Thermostat.ts" />
/// <reference path="Vectors.ts" />
/// <reference path="System.ts" />


namespace Driver {

  export interface IDriver {
    run(s: MolecularMechanics.System): number;
    setCallback(fcn: (arg1: MolecularMechanics.System, arg2: number) => void, callfreq: number): void;
  }
  
  
  export class LeapfrogIntegrator implements IDriver {
    private callbackfreq?: number;
    private callback?: (arg1: MolecularMechanics.System, arg2: number) => void;
    private thermostat: Thermostat.IThermostat;
    
    timestep: number;
    nsteps: number;
    
    constructor(params: {nsteps: number, dt: number}) {
      this.nsteps = params.nsteps;
      this.timestep = params.dt;
    }

    setCallback(fcn: (arg1: MolecularMechanics.System, arg2: number) => void, callfreq: number) {
      this.callback = fcn;
      this.callbackfreq = callfreq;
    }

    setThermostat(thermostat: Thermostat.IThermostat) {
      this.thermostat = thermostat;
    }

    run(system: MolecularMechanics.System): number {
      var masses: number[] = system.getAtoms().map(atom => atom.mass);
      var totalMass = masses.reduce((a, b) => a + b, 0);
      var Ndf = system.countDegreesOfFreedom();
      var velocities = system.getVelocities();
      var coords = system.getCoordinates();
      var forces = system.getForces();
      var Natoms = system.getAtomCount();
      var vtmp = VectorPool.getVector();

      var instantTemp: number;

      var halfdtOverMass: number[] = masses.map(m => 0.5*this.timestep/m);
      // var forces_t: Vector[] = forces.map(v => VectorPool.getVector());
      var forces_t: Vector[] = VectorPool.getVector(Natoms);

      var step = 0;
      var energy = system.eval(coords, forces);
      instantTemp = this.thermostat.getInstantTemperature(velocities,masses,Ndf);
      // console.log(step, energy, instantTemp);
      var vcom = VectorPool.getVector();

      var momRemover = new MomentumRemover({removeLinear: true, removeAngular: true, freq: 10});

      while ((step++) < this.nsteps) {
        
        // temperature and pressure scaling
        this.thermostat.update(velocities, instantTemp);
        
        
        // integrate equations of motion (Verlet velocity)
        for(var i=0; i<Natoms; i++) {
          Vector.omult(forces[i], halfdtOverMass[i], vtmp); // tmp = dt*f_i/(2*m_i)
          vtmp.iadd(velocities[i]);  // tmp += v_i
          vtmp.imult(this.timestep);
          coords[i].iadd(vtmp);
          forces_t[i].set(forces[i].x, forces[i].y, forces[i].z);
        }

        // calculate forces at t+dt
        energy = system.eval(coords, forces);
        // Ekin = this.getKineticEnergy(masses, velocities);
        // instantTemp = 2.0*Ekin/(boltzman*Ndf);
        //instantTemp = this.thermostat.getInstantTemperature(velocities,masses,Ndf);
        
        // update velocities
        for(var i=0; i<Natoms; i++) {
          Vector.oadd(forces_t[i], forces[i], vtmp);
          vtmp.imult(halfdtOverMass[i]);
          velocities[i].iadd(vtmp);
        }

        instantTemp = this.thermostat.getInstantTemperature(velocities,masses,Ndf);
        
        // remove linear/angular momentum
        if (momRemover && (step%10 == 0)) {
          // vcom.izero();
          // velocities.forEach((vel, i) => {
          //   Vector.omult(vel, masses[i], vtmp);
          //   vcom.iadd(vtmp);
          // });
          // vcom.imult(1.0/totalMass);
          // velocities.forEach(v => v.isub(vcom));

          momRemover.apply(coords, velocities, masses);
        }


        if (this.callback && (step%this.callbackfreq == 0)) {
          this.callback(system, step);
          // console.log('>', step, energy, instantTemp);//, velocities[0]);
        }


      }
      VectorPool.collect();
      return -1;
    }
    
  }

  
  
  export class MomentumRemover {
    
    removeAngular: boolean;
    removeLinear: boolean;
    freq: number;

    constructor(prm: {removeAngular: boolean, removeLinear: boolean, freq: number}) {
      this.removeAngular = prm.removeAngular;
      this.removeLinear = prm.removeLinear;
      this.freq = prm.freq;
    }

    apply(xyz: Vector[], vel: Vector[], masses: number[]) {
      var p0: Vector = VectorPool.getVector();    // single particle linear momentum
      var j0: Vector = VectorPool.getVector();    // single particle angular momentum
      var gj: Vector = VectorPool.getVector();    // group angular momentum
      var gx: Vector = VectorPool.getVector();    // group center of mass
      var gv: Vector = VectorPool.getVector();    // group velocity
      var gp: Vector = VectorPool.getVector();    // group linear momentum
      var gw: Vector = VectorPool.getVector();    // group angular velocity
      var dx: Vector = VectorPool.getVector();    // single particle's change in position
      var dv: Vector = VectorPool.getVector();    // single particle's change in velocity
      var gi: Vector[] = VectorPool.getVector(3); // group inertia tensor //[VectorPool.getVector(),VectorPool.getVector(),VectorPool.getVector()];
      var Icm:Vector[] = VectorPool.getVector(3); // group inertia tensor //[VectorPool.getVector(),VectorPool.getVector(),VectorPool.getVector()];
      
      var nAtoms = xyz.length;
      var m0: number;

      for (var i=0; i<nAtoms; i++) {
        m0 = masses[i];
        // linear momentum
        Vector.omult(vel[i], m0, p0);
        gp.iadd(p0);
        
        // angular momentum
        Vector.ocross(xyz[i], vel[i], j0);
        j0.imult(m0);
        gj.iadd(j0);
        gx.set(gx.x + m0*xyz[i].x, gx.y + m0*xyz[i].y, gx.z + m0*xyz[i].z);
        this.update_tensor(xyz[i], m0, gi);
      }
      
      var M = masses.reduce((a, b) => a + b, 0);
      Vector.omult(gp, 1.0/M, gv);

      // calculate group center of mass
      gx.imult(1.0/M);

      // Subtract the center of mass contribution to the angular momentum
      var jcm: Vector = VectorPool.getVector();
      Vector.ocross(gx, gv, jcm);
      jcm.imult(M);
      gj.isub(jcm);

      // Subtract the center of mass contribution from the inertia tensor
      this.update_tensor(gx, M, Icm);
      for(var i=0; i<3; i++)
        gi[i].isub(Icm[i]);

      /* Compute angular velocity, using matrix operation 
      * Since J = I w
      * we have
      * w = I^-1 J
      */
      this.get_minv(gi,Icm);
      gw.set(Icm[0].dot(gj), Icm[1].dot(gj), Icm[2].dot(gj));

      /* Compute the correction to the velocity for each atom */
      for (var i=0; i<nAtoms; i++) {
        vel[i].isub(gv);
        Vector.osub(xyz[i], gx, dx);
        Vector.ocross(gw, dx, dv);
        vel[i].isub(dv);
      }


      p0.free();
      j0.free();
      gj.free();
      gx.free();
      gv.free();
      gp.free();
      gw.free();
      dx.free();
      dv.free();
      while(gi.length>0)
        gi.pop().free();
      while(Icm.length>0)
        Icm.pop().free();
      jcm.free();

    }

    private update_tensor(vec: Vector, s: number, tensor: Vector[]) {
      var xy = vec.x * vec.y * s;
      var xz = vec.x * vec.z * s;
      var yz = vec.y * vec.z * s;
      tensor[0].x += vec.x * vec.x * s;
      tensor[1].y += vec.y * vec.y * s;
      tensor[2].z += vec.z * vec.z * s;
      tensor[0].y += xy;
      tensor[1].x += xy;
      tensor[0].z += xz;
      tensor[2].x += xz;
      tensor[1].z += yz;
      tensor[2].y += yz;
    }

    private get_minv(A: Vector[], B: Vector[]) {
      
      var tmp:Vector[] = VectorPool.getVector(3);

      tmp[0].x =  A[1].y + A[2].z;
      tmp[1].x = -A[0].y;
      tmp[2].x = -A[0].z;
      tmp[0].y = -A[0].y;
      tmp[1].y =  A[0].x + A[2].z;
      tmp[2].y = -A[1].z;
      tmp[0].z = -A[0].z;
      tmp[1].z = -A[1].z;
      tmp[2].z =  A[0].x + A[1].y;

      var rfac: number = (tmp[0].x+tmp[1].y+tmp[2].z)/3;
      var fac = 1/rfac;
      tmp[0].imult(fac);
      tmp[1].imult(fac);
      tmp[2].imult(fac);
      this.m_inv(tmp[0].x,tmp[0].y,tmp[0].z,
                 tmp[1].x,tmp[1].y,tmp[1].z,
                 tmp[2].x,tmp[2].y,tmp[2].z, B);

      B[0].imult(fac);
      B[1].imult(fac);
      B[2].imult(fac);

      while(tmp.length>0)
        tmp.pop().free();
    }


    /**
     * 3x3 matrix inverse
     * @param {number} a{i,j} - 3x3 matrix components
     * @param {Vector[]} inv  - matrix inverse
     */
    private m_inv(a11:number,a12:number,a13:number,
                  a21:number,a22:number,a23:number,
                  a31:number,a32:number,a33:number, inv: Vector[]) {
      var det = a11*(a33*a22-a32*a23) - a21*(a33*a12-a32*a13) + a31*(a23*a12-a22*a13);
      inv[0].set(  a33*a22-a32*a23,  -(a33*a12-a32*a13),   a23*a12-a22*a13);
      inv[1].set(-(a33*a21-a31*a23),   a33*a11-a31*a13,  -(a23*a11-a21*a13));
      inv[2].set(  a32*a21-a31*a22,  -(a32*a11-a31*a12),   a22*a11-a21*a12);
      inv[0].imult(1/det);
      inv[1].imult(1/det);
      inv[2].imult(1/det);
    }

  }
  

  export class SteepestDescent implements IDriver {
    private nSteps: number;
    private Etol: number;
    private stepsize: number;
    private callback?: (arg1: MolecularMechanics.System, arg2: number) => void;

    constructor() {
      this.nSteps = 10000;
      this.Etol = 0.00001;
      this.stepsize = 1;
    }
    setCallback(fcn: (arg1: MolecularMechanics.System, arg2: number) => void) {
      this.callback = fcn;
    }

    run(system: MolecularMechanics.System): number {
      var Eold: number;
      var Enew: number = +Infinity;
      var forces = system.getForces();
      var coords = system.getCoordinates();
      var normSqr: number;
      var stepsize: number;
      
      normSqr = forces.reduce<number>((old, cur, i) => {
        var nSqr = cur.normSq();
        return nSqr>old ? nSqr : old;
      }, -1);

      var step = 0;
      while ((++step) < this.nSteps) {
        // forces.forEach(v => v.izero());
        Eold = Enew;
        Enew = system.eval(coords, forces);

        console.log(step, 'stepsize',this.stepsize, Enew-Eold);
        if (Math.abs(Eold-Enew) < this.Etol)
          break;
        
        if (Enew > Eold)
          this.stepsize *= 0.5
        else
          this.stepsize *= 1.05

        normSqr = forces.reduce<number>((old, cur) => {
          var nSqr = cur.normSq();
          return (nSqr>old ? nSqr : old);
        }, -1);

        stepsize = this.stepsize / Math.sqrt(normSqr);
        
        for (var i=0; i<coords.length; i++) {
          forces[i].imult(stepsize);
          coords[i].iadd(forces[i]);
        }
        if (this.callback)
          this.callback(system, Enew);

      } 
      return Enew;
    }
  }
}