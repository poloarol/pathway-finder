<template>
    <div id="form-job">
        <b-form>
            <b-form-group id='email-input-1' label='Notification Settings' label-for='email-input'>
                <div class="form-row">
                    <div class="col-10">
                        <b-form-input id='email-input' v-model="form.email" type='email' placeholder='E-mail Address' required></b-form-input>
                    </div>
                    <div class="col">
                        <b-icon v-b-tooltip.hover icon='info-circle-fill' title='This allows us to send you the results of the analysis'></b-icon>
                    </div>
                </div>
            </b-form-group>
            <b-form-group label="Data Input">
                <div class="form-row">
                    <div class="col">
                        <div class="btn-toolbar justify-content-between" role="toolbar" aria-label="Toolbar with button groups">
                            <div class="btn-group" role="group" aria-label="First group">
                                <b-button @click="toggle_ncbi_browser(1)" variant='danger' id="upload-file">Upload files</b-button>
                                <b-button @click="toggle_ncbi_browser(2)" variant='outline-danger' id="connect-ncbi">Get from NCBI</b-button>
                            </div>
                        </div>
                    </div>
                    <div class="col-8">
                        <div class="input-group">
                            <div class="container-fluid" id="input_1">
                                <div class="toggler">
                                    <b-form-input v-model="form.accession" placeholder="NCBI GenBank Accession Number"></b-form-input>
                                </div>
                                <b-form-group label-for='file-default' id="input_2">
                                    <b-form-file id="file-default" v-model="form.files" placeholder="Browse" drop-placeholder="Drop file here ..."></b-form-file>
                                </b-form-group>
                            </div>
                        </div>
                    </div>
                </div>
                <!-- <div class="form-row">
                    <div>
                        <b-form-group>
                            <div class="col-4">
                                <b-form-input id='similarity-input' v-model="form.similarity" type='similarity' placeholder='% Similarity' required></b-form-input>
                            </div>
                            <div class="col-4">
                                <b-form-input id='gene-input' v-model="form.gene" type='gene' placeholder='Gene of Interest' required></b-form-input>
                            </div>
                            <div class="col-4">
                                <b-form-input id='basepairs-input' v-model="form.email" type='basepairs' placeholder='Number of Basepairs' required></b-form-input>
                            </div>
                        </b-form-group>
                    </div>
                </div> -->
            </b-form-group>
        </b-form>
    </div>
</template>

<style>
    #form-job {
        width: 80%;
        margin: 0 auto;
    }
    .toggler {
        display: none;
    }
</style>

<script>
import FooterComponentVue from './FooterComponent.vue'
export default {
    name: 'FormComponent',
    data() {
        return {
            form: {}
        }
    },
    methods: {
        toggle_ncbi_browser(number) {
            let upload = document.getElementById('upload-file');
            let ncbi = document.getElementById('connect-ncbi');

            if(number === 1) {
                if(upload.classList.contains('btn-danger'))
                    this.forward_one(upload, ncbi);
                else
                    this.reverse_one(upload, ncbi);
            }else{
                if(ncbi.classList.contains('btn-outline-danger'))
                    this.forward_two(upload, ncbi);
                else
                    this.reverse_two(upload, ncbi);
            }
        },
        forward_one(upload, ncbi){
            upload.classList.remove('btn-danger');
            upload.classList.add('btn-outline-danger');

            ncbi.classList.remove('btn-outline-danger');
            ncbi.classList.add('btn-danger');
        },
        forward_two(upload, ncbi){
            ncbi.classList.remove('btn-outline-danger');
            ncbi.classList.add('btn-danger');

            upload.classList.remove('btn-danger');
            upload.classList.add('btn-outline-danger');
        },
        reverse_one(upload, ncbi){
            upload.classList.remove('btn-outline-danger');
            upload.classList.add('btn-danger');
            
            ncbi.classList.remove('btn-danger');
            ncbi.classList.add('btn-outline-danger');
        },
        reverse_two(upload, ncbi){
            ncbi.classList.remove('btn-danger');
            ncbi.classList.add('btn-outline-danger');
            
            upload.classList.remove('btn-outline-danger');
            upload.classList.add('btn-danger');
        }
    }
}
</script>