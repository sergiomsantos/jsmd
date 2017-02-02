var xmlrpc = require('xmlrpc');

var client = xmlrpc.createClient({ host: 'localhost', port: 9123, path: '/RPC2'})
 
  // client.methodCall('do', ['show sticks'], function (error, value) {
  //   console.log('Method response for \'anAction\': ' + value)
  // });
var pdb = '';
pdb += 'ATOM      1  N   UNK             2.060   3.478   0.676\n';
pdb += 'ATOM      2  C   UNK             3.231   4.154   0.912\n';
pdb += 'ATOM      3  C   UNK             4.110   3.188   1.367\n';
pdb += 'ATOM      4  N   UNK             3.492   1.960   1.409\n';
pdb += 'ATOM      5  C   UNK             2.266   2.166   0.990\n';
pdb += 'ATOM      6  H   UNK             1.196   3.873   0.333\n';
pdb += 'ATOM      7  C   UNK             5.525   3.320   1.781\n';
pdb += 'ATOM      8  H   UNK             5.889   4.344   1.652\n';
pdb += 'ATOM      9  H   UNK             6.158   2.656   1.183\n';
pdb += 'ATOM     10  H   UNK             5.638   3.048   2.835\n';
pdb += 'ATOM     11  C   UNK             3.386   5.607   0.683\n';
pdb += 'ATOM     12  H   UNK             2.447   6.067   0.359\n';
pdb += 'ATOM     13  H   UNK             4.138   5.793  -0.090\n';
pdb += 'ATOM     14  H   UNK             3.707   6.106   1.602\n';
pdb += 'ATOM     15  H   UNK             1.494   1.415   0.893';

//client.methodCall('load', ['/Users/SergioSantos/Dropbox/ccbook/lib/ab.pdb'], function (error, value) {
//  console.log('Method response for \'cmd.load\': ' + value)
//});
client.methodCall('read_pdbstr', [pdb, 'fromnode'], function (error, value) {
  console.log('Method response for \'cmd.read_pdbstr\': ' + value)
});