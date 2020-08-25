<template>
    <div container fluid class='main'>
      <div class='carousel-main'>
        <vue-carousel id='map' v-if='data' :items='1' :nav='true'
          :navigationEnable='true' :center='true'>
          <canvas v-for='(datum, index) in data' :key=index :id="'canvas'+index"
            width='750' height='420'></canvas>
        </vue-carousel>
      </div>
      <div class='carousel-dashboard'>
        <b-card class='card-text'>
          <b-card-text id='org'><b>Organism: </b></b-card-text>
          <b-card-text id='acc'><b>Accession Number: </b></b-card-text>
          <b-card-text id='gene'><b>Gene: </b></b-card-text>
          <b-card-text id='prot'><b>Protein: </b></b-card-text>
          <b-card-text id='loc'><b>Location: </b>Start - Stop</b-card-text>
          <b-button variant="outline-success" class="mb-2" id='protein-file-icon'
            @click="copyToClipBoard(1)">
            <b-icon icon="files" aria-hidden="true"></b-icon> Copy Protein Sequence
          </b-button>
          <b-button variant="outline-success" class="mb-2 push-prot-right" id='protein-download-icon'
            @click="downloadFile(1)">
            <b-icon icon="download" aria-hidden="true"></b-icon> Protein FASTA
          </b-button>
          <b-button variant="outline-info" class="mb-2" id='dna-file-icon'
            @click="copyToClipBoard(2)">
            <b-icon icon="files" aria-hidden="true"></b-icon> Copy DNA Sequence
          </b-button>
          <b-button variant="outline-info" class="mb-2 push-dna-right" id='dna-download-icon'
            @click="downloadFile(2)">
            <b-icon icon="download" aria-hidden="true"></b-icon> DNA FASTA
          </b-button>
          <b-button block variant="outline-dark" class="mb-2">
            <b-icon icon='cloud-download' class="cloud-icon"></b-icon>
            <b-link target='_blank' href="https://www.ncbi.nlm.nih.gov/nuccore/" id="link-modifier">
              <span>Download Pathway</span>
            </b-link>
          </b-button>
        </b-card>
      </div>
    </div>
</template>

<style>
  div.carousel-main {
    width: 800px;
    margin-top: 3%;
    /* margin-right: 10%; */
  }
  div.carousel-dashboard {
    width: 30%;
  }
  .main {
    margin: 0 auto;
    display: flex;
    align-items: center;
    justify-content: center;
  }
  .card-text {
    text-align: left;
  }
  .push-prot-right {
    float: right;
  }
  .push-dna-right {
    float: right;
  }
  a#link-modifier {
    color:black !important;
    text-decoration: none;
  }
  .cloud-icon {
    margin-top: 1% !important;
    margin-right: 0.5em !important;
  }
</style>

<script>
// import jQuery from 'jquery';
import * as Scribl from '../assets/js/scribl.bundle';

export default {
  name: 'Carousel',
  props: ['data'],
  components: {
  },
  data() {
    return {
      slides: [],
      genes: this.data,
      min: 0,
      max: 0,
    };
  },
  methods: {
    copyToClipBoard(val) {
      const divID = val === 1 ? 'protein-file-icon' : 'dna-file-icon';
      const { value } = document.getElementById(divID);
      const input = document.createElement('input');
      input.setAttribute('value', value);
      document.body.appendChild(input);
      input.select();
      const res = document.execCommand('copy');
      document.body.removeChild(input);
      document.execCommand('copy');
      return res;
    },
    downloadFile(val) {
      const divID = val === 1 ? 'protein-file-icon' : 'dna-file-icon';
      const { value } = document.getElementById(divID);
      let accession = document.getElementById('acc').innerHTML;
      accession = accession.split(':');
      const biostring = this.wrap(value, 60);
      const final = `>${accession[1].trim()}\n${biostring}`;
      const fileextname = val === 1 ? 'protein.fasta' : 'dna.fasta';

      const blob = new Blob([final], { type: 'text/csv' });
      if (window.navigator.msSaveOrOpenBlob) {
        window.navigator.msSaveBlob(blob, `${accession[1].trim()}.${fileextname}`);
      } else {
        const elem = window.document.createElement('a');
        elem.href = window.URL.createObjectURL(blob);
        elem.download = `${accession[1].trim()}.${fileextname}`;
        document.body.appendChild(elem);
        elem.click();
        document.body.removeChild(elem);
      }
    },
    wrap(s, w) {
      /*
        eslint max-len: ["error", { "code": 130 }]
      */
      return s.replace(new RegExp(`(?![^\\n]{1,${w}}$)([^\\n]{1,${w}})\\s`, 'g'), '$1\n');
    },
  },
  async updated() {
    this.$nextTick(() => {
      for (let i = 0; i < this.data.length; i += 1) {
        const tmpID = document.querySelector(`#canvas${i}`);
        const tmpCanvas = new Scribl.Scribl(tmpID, 500);
        const track = tmpCanvas.addTrack();
        // let name;
        let genes;
        // let max = 0;
        for (let j = 0; j < this.data[i].length; j += 1) {
          const color = (this.data[i][j].strand === 1) ? '#998ec3' : '#f1a340';
          const { loc } = this.data[i][j];
          const length = Math.abs(loc[1] - loc[0]);
          const { strand } = this.data[i][j];
          const direction = strand === 1 ? '+' : '-';
          const { locus } = this.data[i][j];
          const { gene } = this.data[i][j];
          const { product } = this.data[i][j];
          const protID = this.data[i][j].prot_id;
          const { trans } = this.data[i][j];
          this.min = this.min ? this.min < loc[0] : loc[0];
          this.max = this.max ? this.max > loc[1] : loc[1];
          /*
            eslint max-len: ["error", { "code": 150 }]
          */
          const dna = Object.prototype.hasOwnProperty.call(this.data[i][j], 'dna') ? this.data[i][j].dna : 'No DNA Sequence';
          /*
            eslint max-len: ["error", { "code": 150 }]
          */
          const organism = Object.prototype.hasOwnProperty.call(this.data[i][j], 'org') ? this.data[i][j].dna : ['Homo Sapiens', 'ABCDEF.1'];
          genes = track.addFeature(new Scribl.BlockArrow('track', loc[1], length, direction));
          genes.color = color;
          genes.onMouseover = locus;
          let link = document.getElementById('link-modifier').href;
          link = `${link}${organism[1].trim()}?report=genbank&log$=seqview&from=${this.min}&to=${this.max}`;
          genes.onClick = function hushpuppy() {
            document.getElementById('protein-file-icon').setAttribute('value', trans);
            document.getElementById('dna-file-icon').setAttribute('value', dna);
            document.getElementById('org').innerHTML = `<b>Organism:</b> ${organism[0]} <b>Accession:</b> <em>${organism[1]}</em>`;
            document.getElementById('acc').innerHTML = `<b>Protein Accession Number</b>: ${protID}`;
            document.getElementById('gene').innerHTML = `<b>Gene</b>: ${gene}`;
            document.getElementById('prot').innerHTML = `<b>Protein</b>: ${product}`;
            document.getElementById('loc').innerHTML = `<b>Location</b>: ${loc[0]} - ${loc[1]}`;
            document.getElementById('link-modifier').setAttribute('href', link);
          };
        }
        tmpCanvas.draw();
      }
    });
  },
};
</script>
