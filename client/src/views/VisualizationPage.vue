<template>
  <div>
    <NavBar></NavBar>
    {{ $router.params }}
    <Carousel :data=data></Carousel>
    <Footer></Footer>
  </div>
</template>

<script>
import axios from 'axios';
import NavBar from '@/components/NavBarComponent.vue';
import Footer from '@/components/FooterComponent.vue';
import Carousel from '@/components/CarouselComponent.vue';


export default {
  name: 'Visualization',
  components: {
    NavBar,
    Footer,
    Carousel,
  },
  props: ['value'],
  data() {
    return {
      uuid: this.$router.history.current.params.id,
      data: '',
    };
  },
  methods: {
    getPathway() {
      const main = 'http://localhost:5000/submission';
      const path = `${main}/${this.uuid}`;
      console.log(path);
      axios.post(path, JSON.stringify({ key: this.uuid }))
        .then((res) => {
          this.data = res.data.key;
        })
        .catch((error) => {
          this.data = error;
        });
    },
  },
  created() {
    this.getPathway();
  },
};
</script>
