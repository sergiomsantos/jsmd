/// <reference path="typings/node.d.ts"/>
/// <reference path="System.ts"/>
/// <reference path="Driver.ts"/>

// npm install xmlrpc

var fs = require('fs');
var xmlrpc = require('xmlrpc');


var params = fs.readFileSync('params/ffdat.js', 'utf8');
var client = xmlrpc.createClient({
  host: 'localhost',
  path: '/RPC2',
  port: 9123
});

// var notify = function(data: string) {
//   var fname = 'dat/tmp' + (ii++) +'.xyz';
//   fs.writeFile(fname, data, (err) => {
//     if (err) throw err;
//     client.methodCall('load', [fname, 'molecule'], function (error, value) {
//       // console.log('Method response for \'cmd.load\': ' + value)
//     });
//   });
// }


var fname = 'dat/tmp.xyz';
var notify = function(data: string) {
  fs.appendFileSync(fname, data);
}


var mm = new MolecularMechanics.System();
mm.loadParamsFromJSON(params);

var ab: string[] = fs.readFileSync('params/ab.xyz', 'utf8').split('\n');
//console.log(ab);


var nAtoms = mm.getAtomCount();
var velocities = Vector.getVector(nAtoms);
var forces = Vector.getVector(nAtoms);
var coords = Vector.getVector(nAtoms);

coords.forEach((crd, i) => {
  var xyz = ab[i+2].split(/\s+/);
  crd.set(parseFloat(xyz[1]), parseFloat(xyz[2]), parseFloat(xyz[3]));
});

mm.setForces(forces);
mm.setCoordinates(coords);
mm.setVelocities(velocities);



mm.getAtoms().forEach(atom => console.log(atom));



// var sd = new Driver.SteepestDescent();
// // sd.setCallback(fcn);
// sd.run(mm);
// // mm.eval(mm.getCoordinates(), mm.getForces());

var asXYZ = function(system: MolecularMechanics.System, coords: Vectors.Vector[]): string {
  var n = system.getAtomCount();
  var atoms = system.getAtoms();
  var s = '';
  var lines = coords.map((item, i) => {
    return `${atoms[i].name} ${item.x.toFixed(3)} ${item.y.toFixed(3)} ${item.z.toFixed(3)}`;
  });
  return '  ' + n + '\nTitle\n' + lines.join('\n') + '\n';
}

// notify(asXYZ(mm, coords));
// console.log('number of allocated vectors:',VectorPool.counter);

var fcn = function(system: MolecularMechanics.System, fx: number) {
  //console.log(fx);
  notify(asXYZ(system, system.getCoordinates()));
}


var integrator = new Driver.LeapfrogIntegrator({dt: 0.001, nsteps: 100000});
//var thermostat = new Thermostat.Berendsen({temperature:300, tcoupl: 0.5, dt:0.0001});
var thermostat = new Thermostat.VRescale({temperature:300});
integrator.setCallback(fcn,100);
integrator.setThermostat(thermostat);
thermostat.assignVelocities(mm);
var Ti = thermostat.getInstantTemperature(mm.getVelocities(), mm.getAtoms().map(m=>m.mass), mm.countDegreesOfFreedom());
console.log(Ti);
//console.log(velocities);
integrator.run(mm);


client.methodCall('load', [fname, 'molecule'], function (error, value) {
 console.log('Method response for \'cmd.load\': ' + value)
});


console.log('number of allocated vectors:',VectorPool.counter);