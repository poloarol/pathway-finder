<template>
    <b-container fluid>
        <b-container class="box">
            <b-row col>
                <b-container fluid>
                    <div class="d-flex justify-content-center">
                        <div class="spinner-border text-success" role="status">
                            <span class="sr-only">Loading...</span>
                        </div>
                    </div>
                </b-container>
            </b-row>
            <b-row col-8>
                <b-container fluid class="thank-you">
                    <h5>Thanks for using PathwayFinder</h5>
                    <br>
                    <h5>Your Job can be retrieved using <b><em> {{ key }} </em></b></h5>
                </b-container>
            </b-row>
        </b-container>
    </b-container>
</template>

<style>

    .spinner-border{
        width: 8rem;
        height: 8rem;
        position: absolute;
    }

    .thank-you{
        margin-top: 12em;
    }

    .box{
        margin-top: 5em;
    }

</style>

<script>
export default {
  name: 'ConfirmationMessage',
  props: ['job'],
  data() {
    return {
      key: this.job,
    };
  },
};
</script>
