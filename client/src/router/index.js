import Vue from 'vue';
import VueRouter from 'vue-router';
import Ping from '../components/Ping.vue';

Vue.use(VueRouter);

const router = new VueRouter({
  mode: 'history',
  base: process.env.BASE_URL,
  routes: [
    {
      path: '/ping',
      name: 'Ping',
      component: Ping,
    },
  ],
});

export default router;
