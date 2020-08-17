import Vue from 'vue';
import VueRouter from 'vue-router';
import Pathway from '../views/PathwayFinder.vue';
import Home from '../views/Home.vue';
import Visualization from '../views/VisualizationPage.vue';

Vue.use(VueRouter);

const router = new VueRouter({
  mode: 'history',
  base: process.env.BASE_URL,
  routes: [{
    path: '/pathway',
    name: 'Pathway',
    component: Pathway,
  }, {
    path: '/',
    name: 'Home',
    component: Home,
  }, {
    path: '/submission/:id',
    name: 'Visualization',
    component: Visualization,
  }],
});
export default router;
