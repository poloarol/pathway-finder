<template>
    <div>
        <b-form inline method="POST" id='retrieve-form'
            enctype="multipart/form-data"
             @submit.prevent="onSubmit()">
            <label class="sr-only" for="inline-form-input-name">Name</label>
            <b-input
            id="inline-form-input-name"
            class="mb-2 mr-sm-2 mb-sm-0 col-8"
            placeholder="Job Code"
            ></b-input>
            <b-button type='submit' size='lg' variant="outline-success">Retrieve</b-button>
        </b-form>
    </div>
</template>

<script>
export default {
  name: 'RetrieveJob',
  components: {},
  data() {
    return {
      data: {},
    };
  },
};
</script>
