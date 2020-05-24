<template>
  <div>
    <NavBar></NavBar>
    <div>
      <p>{{ msg }}</p>
    </div>
    <Footer></Footer>
  </div>
</template>

<script>

import axios from 'axios';
import NavBar from '../components/NavBarComponent';
import Footer from '../components/FooterComponent';

export default {
  name: 'Ping',
  components: {
    NavBar,
    Footer
  },
  data(){
    return {
      msg: '',
    };
  },
  methods: {
    getMessage(){
      const path = 'http://localhost:5000/ping';
      axios.get(path)
        .then((res) => {
          this.msg = res.data;
        })
        .catch((error) => {
          console.log(error);
        });
    },
  },
  created() {
    this.getMessage();
  },
};
</script>
