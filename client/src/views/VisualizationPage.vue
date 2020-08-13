<template>
  <div>
    {{ key }}
  </div>
</template>

<script>
import axios from 'axios';


export default {
  name: 'Visualization',
  components: {
  },
  data() {
    return {
      uuid: '',
      data: '',
    };
  },
  methods: {
    created() {
      this.uuid = this.$router.key;
    },
    getPathway() {
      const main = 'http://localhost:5000/diagram';
      const path = `${main}/${this.uuid}`;
      axios.get(path)
        .then((res) => {
          this.data = res.data;
        })
        .catch((error) => {
          this.data = error;
        });
    },
  },
};
</script>
