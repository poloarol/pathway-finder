        for (let j = 0; j < this.data[i].length; j += 1) {
          const color = (this.data[i][j].strand === 1) ? '#998ec3' : '#f1a340';
          const { loc } = this.data[i][j];
          const length = Math.abs(loc[1] - loc[0]);
          const { strand } = this.data[i][j];
          const direction = strand === 1 ? '+' : '-';
          const { locus } = this.data[i][j];
          // const { gene } = this.data[i][j];
          // const { product } = this.data[i][j];
          // const protID = this.data[i][j].prot_id;
          // const { trans } = this.data[i][j];
          this.min = this.min ? this.min < loc[0] : loc[0];
          this.max = this.max ? this.max > loc[1] : loc[1];
          /*
            eslint max-len: ["error", { "code": 150 }]
          */
          // const dna = Object.prototype.hasOwnProperty.call(this.data[i][j], 'dna') ? this.data[i][j].dna : 'No DNA Sequence';
          /*
            eslint max-len: ["error", { "code": 150 }]
          */
          // const organism = Object.prototype.hasOwnProperty.call(this.data[i][j], 'org') ? this.data[i][j].dna : ['Homo Sapiens', 'ABCDEF.1'];
          genes = track.addFeature(new Scribl.BlockArrow('track', loc[1], length, direction));
          genes.color = color;
          genes.onMouseover = locus;
          // let link = document.getElementById('link-modifier').href;
          // link = `${link}${organism[1].trim()}?report=genbank&log$=seqview&from=${this.min}&to=${this.max}`;
          // genes.onClick = function hushpuppy() {
          //   document.getElementById('protein-file-icon').setAttribute('value', trans);
          //   document.getElementById('dna-file-icon').setAttribute('value', dna);
          //   document.getElementById('org').innerHTML = `<b>Organism:</b> ${organism[0]} <b>Accession:</b> <em>${organism[1]}</em>`;
          //   document.getElementById('acc').innerHTML = `<b>Protein Accession Number</b>: ${protID}`;
          //   document.getElementById('gene').innerHTML = `<b>Gene</b>: ${gene}`;
          //   document.getElementById('prot').innerHTML = `<b>Protein</b>: ${product}`;
          //   document.getElementById('loc').innerHTML = `<b>Location</b>: ${loc[0]} - ${loc[1]}`;
          //   document.getElementById('link-modifier').setAttribute('href', link);
          // };
        }