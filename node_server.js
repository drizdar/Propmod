const {spawn} = require('child_process');
const {readFileSync} = require('fs')
const analysis = JSON.parse(readFileSync('data/analysis.json'));
const Propmod = spawn('python', ['sequence/main.py', 'stdin']);
Propmod.stdout.on('data', function(data) {
    console.log(data.toString('utf8'));
});
Propmod.stdin.setEncoding = 'utf-8';
Propmod.stdin.write(JSON.stringify(analysis[0]))
Propmod.stdin.end();
