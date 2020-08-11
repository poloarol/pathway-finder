<template>
  <div>
    <NavBar></NavBar>
    <div>
      <Form class="show-on-submit" id="form-submit" @gotoggle='toggle'></Form>
      <Confirmation msg="" job="" class="hide-on-submit" id="message-submit"></Confirmation>
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
      msg: '',
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
    toggle() {
      const form = document.getElementById('form-submit');
      const msg = document.getElementById('message-submit');
      form.classList.remove('show-on-submit');
      form.classList.add('hide-on-submit');
      msg.classList.remove('hide-on-submit');
      msg.classList.add('show-on-submit');
    },
  },
  created() {
    this.getMessage();
  },
};
</script>
