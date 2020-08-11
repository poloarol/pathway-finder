/**
 * Genbank parser
 */
function genbank(file, bchart) {
  const lines = file.split('\n');
  const re = new RegExp(/\s+gene\s+([a-z]*)\(?(\d+)\.\.(\d+)/);
  const genes = [];
  let max;
  let min;

  // parse genes
  for (let j = 0; j < lines.length; j++) {
    var gene_info;
    if (gene_info = lines[j].match(re)) {
      gene_info.shift();
      genes.push(gene_info);

      // determine scale dimensions
      var end = gene_info[2];
      if (max == undefined || max > end) max = end;
      var position = gene_info[1];
      if (min == undefined || min < position) min = position;
    }
  }


  // set scale dimensions
  bchart.scale.max = max;
  bchart.scale.min = min;

  // add genes to chart
  for (let i = 0; i < genes.length; i++) {
    // get positional values
    let strand = '+';
    if (genes[i][0] == 'complement') strand = '-';

    var position = genes[i][1];
    var end = genes[i][2];
    position = position - 1 + 1; // force to be integer - TODO make bChart catch non-ints automatically and gracefully fail
    const length = end - position;

    bchart.addGene(position, length, strand);
  }
}
