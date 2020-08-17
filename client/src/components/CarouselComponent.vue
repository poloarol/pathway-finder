<template>
    <div container fluid>
      <!-- <vue-carousel id='map' v-if="data"
        :center="true" :items="1" :nav="true" :navigationEnable="true">
        <div v-for="(datum, index) in data" :key="index"> {{ datum }}</div>
      </vue-carousel> -->
      <ul v-if='data'>
        <li v-for='(datum, index) in data[0]' :key='index'> {{ datum }}</li>
      </ul>
    </div>
</template>

<style>
  canvas {
    margin: 0 auto;
    padding: 10em;
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
    };
  },
  methods: {
    generateDiagram(chart, data) {
      let gene;
      for (let i = 0; i < 1; i += 1) {
        const track = chart.addTrack();
        let name;
        for (let j = 0; j < data[i].length; j += 1) {
          const color = (data[i][j].strand === '+') ? '#998ec3' : '#f1a340';
          gene = track.addFeature(new Scribl.BlockArrow('complex', data[i][j].start, data[i][j].length, data[i][j].strand, color));
          const tmp = data[i][j].extraclass[1] !== '' ? data[i][j].extraclass[1] : data[i][j].extraclass[0];
          name = data[i][j].name !== '' ? data[i][j].name : tmp;
          gene.onMouseover = name;
        }
      }
      chart.draw();
    },
  },
  async updated() {
    // this.$nextTick(() => {
    //   for (let i = 0; i < this.data.length; i += 1) {
    //     const tmpID = document.querySelector(`#canvas${i}`);
    //     const tmpCanvas = new Scribl.Scribl(tmpID, 500);
    //     const track = tmpCanvas.addTrack();
    //     // let name;
    //     let genes;
    //     let max = 0;
    //     for (let j = 0; j < this.data[i].length; j += 1) {
    //       const color = (this.data[i][j].strand === 1) ? '#998ec3' : '#f1a340';
    //       const { loc } = this.data[i][j];
    //       const length = Math.abs(loc[1] - loc[0]);
    //       const { strand } = this.data[i][j];
    //       // const { locus } = this.data[i][j];
    //       // const { gene } = this.data[i][j];
    //       // const { product } = this.data[i][j];
    //       // const protID = this.data[i][j].prot_id;
    //       max = length + 500;
    //       if (j === 0) {
    //         genes = track.addFeature(new Scribl.BlockArrow('track', loc[1], length, strand));
    //       } else {
    //         genes = track.addFeature(new Scribl.BlockArrow('track', max, length, strand));
    //       }
    //       genes.color = color;
    //     }
    //     track.draw();
    //   }
    // });
  },
};
</script>
