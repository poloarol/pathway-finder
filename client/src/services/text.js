export default {
  getNavBarTitle() {
    return [{
      value: 'Home',
      icon: 'shop',
      link: '/',
    },
    {
      value: 'Download - Pathway Finder',
      icon: 'cloud-download',
      link: '/download',
      target: '_target',
    },
    {
      value: 'BGC-Finder',
      icon: 'bounding-box',
      link: '/bgc-finder',
      target: '_target',
    },
    {
      value: 'About',
      icon: 'info-circle-fill',
      link: '/about',
      target: '_target',
    },
    ];
  },
  getNavLinks() {
    return [{
      icon: 'github',
      link: 'https://github.com/poloarol/pathway-finder',
    },
    {
      icon: 'diamond',
      link: 'http://www.boddylab.ca/',
    },
    {
      icon: 'diamond',
      link: 'https://twitter.com/boddylab?lang=en',
    },
    ];
  },
};
