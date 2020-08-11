<template>
  <b-container fluid>
      <NavBar/>
      <Footer/>
  </b-container>
</template>

<script>
// @ is an alias to /src
import NavBar from '@/components/NavBarComponent.vue';
import Footer from '@/components/FooterComponent.vue';

export default {
  name: 'Home',
  components: {
    NavBar,
    Footer,
  },
};
</script>
