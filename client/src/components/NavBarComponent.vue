<template>
    <div>
        <b-navbar role='navigation' id="navbar-whole">
            <b-navbar-brand class="brand-nav" tag="router-link" :to="{path: '/'}" target="_target">
                <img src="../assets/penicillin_chemdraw_one.png" width="150"
                    class="d-inline-block align-top" alt="Penicillin Molecule">
            </b-navbar-brand>

            <b-navbar-toggle target="nav-collapse"></b-navbar-toggle>

            <b-collapse is-nav>
                <b-navbar-nav right id="navigation">
                    <b-nav-item v-for="(title, index) in titles"
                        :key=index :href=title.link id="nav-item"
                        :target=title.target>
                        <b-icon :icon=title.icon font-scale="1.5"></b-icon><br>
                        <span>{{ title.value }}</span>
                    </b-nav-item>
                </b-navbar-nav>
            </b-collapse>

        </b-navbar>
    </div>
</template>

<style>
    #navbar-whole {
        background-color: #8F001A;
        min-height: 4rem;
        padding: 1.5em;
    }
    a{
        color: white !important;
    }
    #nav-item{
        padding: 1.5em !important;
    }
    #navigation{
        margin-left: 33.5em;
    }
    .brand-nav{
        margin-left: 5rem;
    }
</style>

<script>

import texts from '../services/text';

export default {
  name: 'NavBarComponent',
  data() {
    return {
      titles: texts.getNavBarTitle(),
      icons: texts.getNavLinks(),
    };
  },
};
</script>
