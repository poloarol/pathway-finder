<template>
    <nav class="navbar fixed-bottom navbar-dark bg-dark navbar-expand-sm">
        <div class="container-fluid">
            <b-navbar-nav center id="navbar-footer">
                <div class="container-fluid" id="img-logo">
                    <!-- <b-img src="../assets/uottawa_ver.png" fluid alt="uOttawa Logo" width="200" height="auto"></b-img> -->
                    <img src="../assets/uottawa_ver.png" class="img-fluid" alt="uOttawa logo" height="auto" width="80">
                </div>
                <b-nav-item v-for="(icon, index) in icons" :key=index :href="`${icon.link}`" target='_blank' id=footer-nav-item>
                    <b-icon :icon=icon.icon></b-icon>
                </b-nav-item>
            </b-navbar-nav>
        </div>
    </nav>
</template>

<style>

    #navbar-footer{
        margin: 0 auto;
    }
    #footer-nav-item{
        padding: 2em;
    }
    #img-logo{
        margin: 0 auto;
    }
</style>

<script>

import texts from '../services/text'

export default {
    name: 'FooterComponent',
    data() {
        return {
            icons: texts.getNavLinks()
        };
    }
    
}
</script>