<template>
  <ConfirmationMessage :job=uuid />
</template>

<style>

</style>

<script>

import ConfirmationMessage from '@/components/ConfirmationMessage.vue';

export default {
  name: 'Confirmation',
  props: ['job'],
  components: {
    ConfirmationMessage,
  },
  data() {
    return {
      uuid: this.job,
    };
  },
};
</script>
