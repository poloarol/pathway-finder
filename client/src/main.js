import Vue from 'vue';
import { BootstrapVue, IconsPlugin, BootstrapVueIcons } from 'bootstrap-vue';
import { library } from '@fortawesome/fontawesome-svg-core';
import { faGithub, faTwitter } from '@fortawesome/free-brands-svg-icons';
import { faSpinner } from '@fortawesome/free-solid-svg-icons';
// import { far } from '@fortawesome/free-regular-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/vue-fontawesome';
import VueCarousel from 'vue-owl-carousel';


import 'bootstrap/dist/css/bootstrap.css';
import 'bootstrap-vue/dist/bootstrap-vue.css';


import App from './App.vue';
import router from './router';

library.add(faGithub, faTwitter, faSpinner);

// Install BootstrapVue
Vue.use(BootstrapVue);
Vue.use(BootstrapVueIcons);
// Optionally install the BootstrapVue icon components plugin
Vue.use(IconsPlugin);

Vue.component('font-awesome-icon', FontAwesomeIcon);
Vue.component('vue-carousel', VueCarousel);

Vue.config.productionTip = false;

new Vue({
  router,
  render: (h) => h(App),
}).$mount('#app');
