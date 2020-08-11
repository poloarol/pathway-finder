/**
 * Bed parser
 */
function bed(file, chart) {
  const lines = file.split('\n');
  const features = [];
  const max;
  const min;

  // 	track name=pairedReads description="Clone Paired Reads" useScore=1
  const trackInfo = lines[0];

  // parse genes
  numFeatures = lines.length;
  for (let j = 1; j < numFeatures; j++) {
    if (lines[j] == '') break;

    const fields = lines[j].split(' ');

    // chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
    const chromStart = parseInt(fields[1]);
    const chromEnd = parseInt(fields[2]);
    const name = `${fields[0]}: ${fields[3]}`;
    const orientation = fields[5];
    const itemRgb = fields[8];
    const blockLengths = fields[10].split(',');
    const blockStarts = fields[11].split(',');

    const complex = chart.addFeature(new Complex('complex', chromStart, chromEnd, orientation, [], { color: itemRgb, name }));

    for (let k = 0; k < blockLengths.length; k++) {
      if (blockLengths[k] == '') break;
      complex.addSubFeature(new BlockArrow('complex', parseInt(blockStarts[k]), parseInt(blockLengths[k]), orientation));
    }
  }
}
