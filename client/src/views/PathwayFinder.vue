<template>
  <div>
    <NavBar></NavBar>
    <div>
      <Form id="form-submit" v-if="!uuid" @gotoggle='toggle'></Form>
      <Confirmation :job=uuid v-if="uuid" id="message-submit"></Confirmation>
    </div>
    <Footer></Footer>
  </div>
</template>

<style>

  #form-submit{
    margin-top: 5em;
    margin-left: 10em;
  }

  .hide-on-submit{
    display: none;
  }

  .show-on-submit{
    display: block;
  }

</style>

<script>

import axios from 'axios';
import NavBar from '@/components/NavBarComponent.vue';
import Footer from '@/components/FooterComponent.vue';
import Form from '@/components/FormComponent.vue';
import Confirmation from './ConfirmationPage.vue';

export default {
  name: 'Pathway',
  components: {
    NavBar,
    Footer,
    Form,
    Confirmation,
  },
  data() {
    return {
      uuid: '',
    };
  },
  methods: {
    getMessage() {
      const path = 'http://localhost:5000/ping';
      axios.get(path)
        .then((res) => {
          this.msg = res.data;
        })
        .catch((error) => {
          this.msg = error;
        });
    },
    toggle(value) {
      // const form = document.getElementById('form-submit');
      // const msg = document.getElementById('message-submit');
      // form.classList.remove('show-on-submit');
      // form.classList.add('hide-on-submit');
      // msg.classList.remove('hide-on-submit');
      // msg.classList.add('show-on-submit');
      this.uuid = value;
    },
  },
  created() {
    this.getMessage();
  },
};
</script>
