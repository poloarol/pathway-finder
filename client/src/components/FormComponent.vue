<template>
  <div id="form-job">
    <b-form id="form-NCBI" method="POST"
      enctype="multipart/form-data"
      @submit.prevent="onSubmit()">
        <b-form-group id='email-input-1'
          label='Notification Settings'
          label-for='email-input'>
          <div class="form-row">
            <div class="col-10">
              <b-form-input
                id='email-input'
                v-model="form.email"
                type='email'
                placeholder='E-mail Address'>
              </b-form-input>
            </div>
            <div class="col">
              <b-icon v-b-tooltip.hover
                icon='info-circle-fill'
                title='This allows us to send you the results of the analysis'>
              </b-icon>
            </div>
         </div>
        </b-form-group>
        <b-form-group label="BGC Information - From NCBI">
          <div class="form-row">
            <div class="col-6">
              <b-form-input
                id='accession-input'
                type="text"
                v-model="form.accession"
                placeholder='Accession number'>
              </b-form-input>
            </div>
            <div class="col-6">
              <b-form-input
                id='protein-input'
                type="text"
                v-model="form.protein"
                placeholder='Core protein'>
              </b-form-input>
            </div>
          </div>
        </b-form-group>
        <b-form-group label='BGC information - From Sequence'>
          <div class='form-row'>
            <div class='col-12'>
              <b-form-textarea
                id='sequence-input'
                type='text'
                v-model='form.sequence'
                placeholder='Protein Sequence'
                rows='4'>
              </b-form-textarea>
            </div>
          </div>
        </b-form-group>
        <b-form-group label="Detection Parameters">
          <div class="form-row">
            <div class="col-6">
              <b-form-input
                id='similarity-input'
                type="number"
                v-model="form.similarity"
                placeholder='% Similarity'>
              </b-form-input>
            </div>
            <div class="col-6">
              <b-form-input
                id='basepairs-input'
                type="number"
                v-model="form.basepairs"
                placeholder='Number of Basepairs'>
              </b-form-input>
            </div>
            <!-- <div class="col-4">
              <b-form-input
                id='e-value-input'
                type='number'
                v-model='form.evalue'
                placeholder='10'>
              </b-form-input>
            </div> -->
          </div>
        </b-form-group>
      <br><br>
      <b-button type="submit" size='lg' variant="outline-success">Submit</b-button>
    </b-form>
  </div>
</template>

<style>
    #form-job {
        width: 100%;
        display: flex;
        justify-content: center;
        align-items: center;
    }
    #form-NCBI{
        width: 90%;
        text-align: left;
    }
    .toggler {
        display: none;
    }
</style>

<script>

import axios from 'axios';
import { v4 as uuid4 } from 'uuid';

export default {
  name: 'FormComponent',
  data() {
    return {
      form: {},
    };
  },
  methods: {
    toggle_NCBI_browser(number) {
      const UPLOAD = document.getElementById('upload-file');
      const NCBI = document.getElementById('connect-ncbi');

      const UPLOAD_INPUT = document.getElementById('input_1');
      const NCBI_INPUT = document.getElementById('input_2');

      if (number === 1) {
        if (UPLOAD.classList.contains('btn-danger')) {
          this.forward_one(UPLOAD, NCBI, UPLOAD_INPUT, NCBI_INPUT);
        } else {
          this.reverse_one(UPLOAD, NCBI, UPLOAD_INPUT, NCBI_INPUT);
        }
      } else {
        if (NCBI.classList.contains('btn-outline-danger')) {
          this.forward_two(UPLOAD, NCBI, UPLOAD_INPUT, NCBI_INPUT);
        } else {
          this.reverse_two(UPLOAD, NCBI, UPLOAD_INPUT, NCBI_INPUT);
        }
        this.doSomething();
      }
    },
    forward_one(upload, ncbi, up1, up2) {
      upload.classList.remove('btn-danger');
      upload.classList.add('btn-outline-danger');
      ncbi.classList.remove('btn-outline-danger');
      ncbi.classList.add('btn-danger');
      up1.classList.remove('toggler');
      up2.classList.add('toggler');
    },
    forward_two(upload, ncbi, up1, up2) {
      ncbi.classList.remove('btn-outline-danger');
      ncbi.classList.add('btn-danger');
      upload.classList.remove('btn-danger');
      upload.classList.add('btn-outline-danger');
      up1.classList.remove('toggler');
      up2.classList.add('toggler');
    },
    reverse_one(upload, ncbi, up1, up2) {
      upload.classList.remove('btn-outline-danger');
      upload.classList.add('btn-danger');
      ncbi.classList.remove('btn-danger');
      ncbi.classList.add('btn-outline-danger');
      up2.classList.remove('toggler');
      up1.classList.add('toggler');
    },
    reverse_two(upload, ncbi, up1, up2) {
      ncbi.classList.remove('btn-danger');
      ncbi.classList.add('btn-outline-danger');
      upload.classList.remove('btn-outline-danger');
      upload.classList.add('btn-danger');
      up2.classList.remove('toggler');
      up1.classList.add('toggler');
    },
    doSomething() {
      return false;
    },
    msg_form_toggler(key) {
      this.$emit('gotoggle', key);
    },
    onSubmit() {
      const path = 'http://localhost:5000';
      this.form.uuid4 = uuid4();
      this.msg_form_toggler(this.form.uuid4);
      axios.post(`${path}/pathway`, JSON.stringify(this.form))
        .then((res) => {
          if (res.status === 200) {
            this.$router.push({ path: `/submission/${this.form.uuid4}`, params: { value: this.form.uuid4 } });
          }
        });
    },
  },
};
</script>
