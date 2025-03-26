
import React from 'react';
import { Link } from 'react-router-dom';

const Footer: React.FC = () => {
  return (
    <footer className="p-4 bg-gradient-to-r from-gray-800 to-black text-white flex items-center justify-between animate-fadeIn" style={{ animationDelay: '500ms' }}>
      <div>
        <p className="text-gray-300 text-sm">&copy; 2025 DrugDiscov. All rights reserved.</p>
      </div>
      <div className="flex space-x-6">
        <Link to="/privacy" className="text-gray-300 hover:text-blue-400 text-sm">
          Privacy Policy
        </Link>
        <Link to="/terms" className="text-gray-300 hover:text-blue-400 text-sm">
          Terms of Service
        </Link>
        <Link to="/contact" className="text-gray-300 hover:text-blue-400 text-sm">
          Contact Us
        </Link>
      </div>
    </footer>
  );
};

export default Footer;
