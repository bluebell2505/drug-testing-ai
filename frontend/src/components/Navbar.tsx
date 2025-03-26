
import React from 'react';
import { Link } from 'react-router-dom';
import { cn } from '@/lib/utils';
import { useToast } from "@/hooks/use-toast";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuLabel,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu";

const Navbar: React.FC = () => {
  const { toast } = useToast();
  
  // Mock user data
  const userName = "Dr. Smith";
  const userInitials = "DS";

  const handleLogout = () => {
    console.log('User logged out');
    toast({
      title: "Success",
      description: "Logged out successfully!",
    });
    // In a real app, would redirect to login page
  };

  return (
    <nav className="flex items-center justify-between px-6 py-4 bg-gradient-to-r from-gray-800 to-black">
      <Link to="/dashboard" className="flex items-center gap-1 text-2xl font-bold">
        <span className="text-blue-400">Drug</span>
        <span className="text-white">Discov</span>
      </Link>
      
      <div className="hidden lg:flex items-center space-x-8">
        <Link to="/dashboard" className="text-white hover:text-blue-400">
          Dashboard
        </Link>
        <Link to="/projects" className="text-white hover:text-blue-400">
          Projects
        </Link>
        <Link to="/molecule-generator" className="text-white hover:text-blue-400">
          Molecule Generator
        </Link>
        <Link to="/analytics" className="text-white hover:text-blue-400">
          Analytics
        </Link>
      </div>
      
      <DropdownMenu>
        <DropdownMenuTrigger className="flex items-center gap-2 outline-none">
          <div className="bg-blue-500 text-white w-10 h-10 flex items-center justify-center rounded-full font-medium">
            {userInitials}
          </div>
        </DropdownMenuTrigger>
        <DropdownMenuContent align="end" className="w-56 bg-gray-800 border border-gray-700 text-white">
          <DropdownMenuLabel>My Account</DropdownMenuLabel>
          <DropdownMenuSeparator className="bg-gray-700" />
          <DropdownMenuItem className="hover:bg-gray-700 cursor-pointer">
            <Link to="/profile" className="flex items-center gap-2 w-full">
              Profile
            </Link>
          </DropdownMenuItem>
          <DropdownMenuItem className="hover:bg-gray-700 cursor-pointer">
            <Link to="/settings" className="flex items-center gap-2 w-full">
              Settings
            </Link>
          </DropdownMenuItem>
          <DropdownMenuSeparator className="bg-gray-700" />
          <DropdownMenuItem 
            className="hover:bg-gray-700 text-red-400 cursor-pointer"
            onClick={handleLogout}
          >
            <div className="flex items-center gap-2 w-full">
              Logout
            </div>
          </DropdownMenuItem>
        </DropdownMenuContent>
      </DropdownMenu>
    </nav>
  );
};

export default Navbar;
