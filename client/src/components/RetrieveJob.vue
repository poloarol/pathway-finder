<template>
  <div>
    <b-form inline method="POST"
      enctype="multipart/form-data"
      @submit.prevent="onSubmit()">
      <label class="sr-only" for="inline-form-input-name">Name</label>
      <b-input
        id="inline-form-input-name"
        class="mb-2 mr-sm-2 mb-sm-0 col-10"
        placeholder="Job Code"
        v-model='form.jobcode'>
      </b-input>
      <b-button type='submit' size='lg' variant="outline-success">Retrieve</b-button>
    </b-form>
  </div>
</template>


<script>
import axios from 'axios';

export default {
  name: 'RetrieveJob',
  components: {},
  data() {
    return {
      form: {},
    };
  },
  methods: {
    onSubmit() {
      const path = 'http://localhost:5000/submission';
      axios.post(`${path}/${this.form.jobcode}`, JSON.stringify({ key: this.form.jobcode }))
        .then((res) => {
          if (res.status === 200) {
            this.$router.push({ path: `/submission/${this.form.jobcode}`, params: { value: this.form.jobcode } });
          }
        });
    },
  },
};
</script>
